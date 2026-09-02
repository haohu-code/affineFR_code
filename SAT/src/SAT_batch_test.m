% SAT_BATCH_TEST  Batch test script for SAT problem relaxations.
%
%   [EXPoutputs, folderNames] = SAT_batch_test(options) runs a series of SAT
%   problem tests over all datasets found in the SAT data directory.
%
%   Inputs:
%       options - Struct containing configuration flags:
%           .mode       - Execution mode (1: FRA only, 2+: Full relaxations)
%           .obj_mode   - Objective type ('linear', 'quadratic')
%           .robust     - Boolean; if true, runs SDP in a separate MATLAB process
%           .relaxation - Type of SDP relaxation to use
%           .computeV   - Boolean; construct V during FRA (default true)
%           .datasets   - Optional cell array of dataset folder names
%           .instances  - Optional cell array of instance names to run
%           .useTable6Objectives - Load the fixed Table 6 linear objectives
%           .skipOriginalSDPInstances - Instances not to run for safety
%
%   Outputs:
%       EXPoutputs - Cell array containing results for each dataset.
%       folderNames - Dataset names in the same order as EXPoutputs.

function [EXPoutputs,folderNames] = SAT_batch_test(options)


% Get data folder names.
if isfield(options,'datasets')
    folderNames = cellstr(options.datasets(:));
else
    folderNames = get_data_folder();
end
m = length(folderNames);

% initialize the output
EXPoutputs = cell(m,1);
for i = 1:m
    fprintf('\n\n\nThe %d-th dataset : %s \n',i,folderNames{i})
    EXPoutputs{i} = SAT_test(folderNames{i},options);
    fprintf('\n\n')
end



end



% SAT_TEST  Runs relaxations on a specific SAT dataset.
%
%   outputs = SAT_test(data_set, options) processes all .cnf files in
%   the specified data_set folder.

function outputs = SAT_test(data_set,options)

% load the options
mode = options.mode;
if isfield(options,'obj_mode')
    obj_mode = options.obj_mode;
else
    obj_mode = 'none';
end
robust = options.robust;
if isfield(options,'computeV')
    computeV = options.computeV;
else
    computeV = true;
end


if robust
    [cmd,cmd_fr] = get_cmd();
end

% --- Solver Size Limits ---
% Skip instances larger than these thresholds
MAX_SIZE = [100000;150000]; % [variables; clauses]

% Maximum size for SDP solvers to avoid memory overflow or hung processes
MAX_SDP_SIZE = [350;2000];

% Maximum size for Integer Linear Programming (ILP) solvers
MAX_ILP_SIZE = [1500;4000];


% load the list of problem instances in the data_set
fileList = get_files(data_set);
if isfield(options,'instances')
    requestedNames = cellstr(string(options.instances(:)));
    fileNames = regexprep({fileList.name},'\.cnf$','');
    fileList = fileList(ismember(fileNames,requestedNames));
end
K = length(fileList);

% print the title
print_title(mode)

% loop through the problem
outputs = cell(K,1);
for i = 1:K
    % Correct way to join the folder and name from the dir struct
    currentFile = fullfile(fileList(i).folder, fileList(i).name);
    % name = fileList(i).name;
    % if ~strcmp(name(1:end-4),'jnh210')
    %     continue
    % end

    % keyboard

    % Load the sat instance
    satProblem = loadSATLIB(currentFile);

    % check the problem size
    n = satProblem.numVars;
    % m = satProblem.numClauses; % (this is different from below.)
    m = length(satProblem.clauses);

    % Skip instances that are too large for any processing
    if or(n > MAX_SIZE(1),m > MAX_SIZE(2))
        continue
    end

    % Skip instances too large for SDP if in SDP mode
    if and(mode >= 2,or(n > 5*MAX_SDP_SIZE(1),m > MAX_SDP_SIZE(2)))
        continue
    end

    % construct the ILP formulation
    model = SAT_ILP_FEA(satProblem,obj_mode);
    if isfield(options,'useTable6Objectives') && ...
            options.useTable6Objectives
        [fixedObjective,objectiveRecord] = ...
            load_table6_objective(fileList(i).name);
        if objectiveRecord.numVariables ~= n
            error('Stored objective length does not match %s.', ...
                fileList(i).name);
        end
        model.obj = fixedObjective;
    end

    % save the instance information
    outputs{i}.info = save_instance(fileList(i).name,model,m,n);

    % print the instance information
    print_instance(outputs{i}.info)

    % --- Affine Facial Reduction (FRA) ---
    % Apply affine facial reduction analysis to reduce problem dimension
    outputs{i}.affineFRA = FRAcompute(model,'affineFR',computeV);
    print_FRA(outputs{i}.affineFRA)
    r = outputs{i}.affineFRA.r; % Reduced dimension
    N = outputs{i}.affineFRA.N;

    % remove some memory consuming variables
    outputs{i}.affineFRA = rmfield(outputs{i}.affineFRA,'W');
    outputs{i}.affineFRA = rmfield(outputs{i}.affineFRA,'V');
    outputs{i}.affineFRA = rmfield(outputs{i}.affineFRA,'N');
    outputs{i}.affineFRA = rmfield(outputs{i}.affineFRA,'solver_output');

    % If mode is 1, we only care about FRA results
    if mode == 1
        fprintf('\n')
        if isfield(outputs{i}.info,'Q')
            outputs{i}.info = rmfield(outputs{i}.info,'Q');
        end
        outputs{i}.info = rmfield(outputs{i}.info,'obj');
        continue
    end

    % Skip further computations if FRA achieved negligible reduction (too dense)
    if r >= (n+1)*.95
        fprintf('\n')
        continue
    end

    % if the FR reduced instance is still too big, then skip it.
    if or(r > MAX_SDP_SIZE(1),m > MAX_SDP_SIZE(2))
        fprintf('\n')
        continue
    end

    % --- Solve ILP ---
    params.outputflag = 0; % Set to 0 to disable solver logs
    params.TimeLimit = 600; % Sets limit to 10 minutes
    outputs{i}.ILP = [];
    if and(n <= MAX_ILP_SIZE(1),m <= MAX_ILP_SIZE(2))
        outputs{i}.ILP = gurobi(model, params);
    end
    print_result(outputs{i}.ILP)

    if strcmp(obj_mode,'linear')
        % --- Solve LP relaxation ---
        model.vtype = 'C';
        outputs{i}.LP = gurobi(model, params);
        print_result(outputs{i}.LP)
    elseif strcmp(obj_mode,'quadratic')
        % solve the RLT relaxation
        outputs{i}.RLT = [];
        if all([r <= MAX_SDP_SIZE(1),m <= MAX_SDP_SIZE(2)])
            output = SAT_relaxation(model,'RLT');
            outputs{i}.RLT = output;
        end
        print_result(outputs{i}.RLT)
    end


    % --- SDP Relaxation ---
    % Solve standard Semidefinite Programming (SDP) relaxation
    outputs{i}.SDP = [];
    skipOriginalSDP = false;
    if isfield(options,'skipOriginalSDPInstances')
        skipNames = cellstr(string(options.skipOriginalSDPInstances(:)));
        skipOriginalSDP = ismember(outputs{i}.info.name,skipNames);
    end
    if skipOriginalSDP
        outputs{i}.SDP = struct( ...
            'status','SKIPPED_OOM_SAFETY', ...
            'historicalStatus','OOM', ...
            'objval',[], ...
            'runtime',[]);
    elseif all([n+1 <= MAX_SDP_SIZE(1),m <= MAX_SDP_SIZE(2)])
        relaxation = options.relaxation;
        if robust == true
            % Use external MATLAB call for robustness against crashes
            save('robust_sdp_run.mat','model','relaxation')
            system(cmd);
            load('robust_sdp_run.mat','output')
            outputs{i}.SDP = output;
        else
            outputs{i}.SDP = SAT_relaxation(model,relaxation);
        end
    end
    print_result(outputs{i}.SDP)


    % --- SDP Facial Reduction (SDPFR) ---
    % Solve facial-reduced SDP relaxation (lower dimension)
    outputs{i}.SDPFR = [];
    if all([r <= MAX_SDP_SIZE(1),m <= MAX_SDP_SIZE(2)])
        relaxation = options.relaxation;
        if robust == true
            % Use external MATLAB call for facial-reduced version
            save('robust_sdp_run.mat','model','relaxation','N')
            system(cmd_fr);
            load('robust_sdp_run.mat','output')
            outputs{i}.SDPFR = output;
        else
            outputs{i}.SDPFR = SAT_FR_relaxation(model,1,N);
        end
    end
    print_result(outputs{i}.SDPFR)

    fprintf('\n')
end

end

function fileList = get_files(data_set)
% GET_FILES  Retrieves all .cnf files within a specific dataset folder.

projectRoot = get_project_root();
dataDir = fullfile(projectRoot, 'data', data_set);

% Now look for the files
fileList = dir(fullfile(dataDir,'**', '*.cnf'));

% Diagnostic: If still empty, check the path string
if isempty(fileList)
    fprintf('Path being searched: %s\n', dataDir);
    error('No files found! Check if the folder name is correct.');
end

end

function print_title(mode)
% PRINT_TITLE  Displays the table header for the test results.


fprintf('%-26s%8s%8s','Name','Var','Clause')
fprintf('%3s',' | ')
fprintf('%8s%8s%16s','r','time','status')
fprintf('%3s',' | ')
if mode >= 2
    fprintf('%8s%8s%16s','ILP_obj','time','status')
    fprintf('%3s',' | ')
    fprintf('%8s%8s%16s','LP_obj','time','status')
    fprintf('%3s',' | ')
    fprintf('%8s%8s%16s%8s%8s','SDP_obj','time','status')
    fprintf('%3s',' | ')
    fprintf('%8s%8s%16s%8s%8s','SDPFR_obj','time','status','vtime')
    fprintf('%3s',' | ')
end
fprintf('\n')

end



function print_instance(output)
% print the instance name, size
fprintf('%-26s%8d%8d',output.name,output.n,output.m)
fprintf('%3s',' | ')
end



function print_FRA(output)
% prin the affine FRA results
fprintf('%8d%8.2f%16s',output.r,output.total_time,output.status)
fprintf('%3s',' | ')
end


function print_result(output)
% PRINT_RESULT  Generic formatter for solver results (objective, time, status).

% Handle missing or incomplete output structures
if ~isfield(output,'status')
    output.status = [];
end

if strcmp(output.status,'INFEASIBLE')
    output.objval = -inf;
end

if ~isfield(output,'objval')
    output.objval = [];
end

if ~isfield(output,'runtime')
    output.runtime = [];
end

if ~isfield(output,'mu')
    output.mu = -1;
end

if ~isfield(output,'pfeas_total')
    output.pfeas_total = [];
end


% Print values with fixed-width formatting for the table
fprintf('%8.2f%8.2f%16s',output.objval,output.runtime,output.status)
if isfield(output,'vtime')
    fprintf('%8.2f',output.vtime)
end
fprintf('%3s',' | ')

end


function folderNames = get_data_folder()
% GET_DATA_FOLDER  Scans the SAT data directory for dataset subfolders.

projectRoot = get_project_root();
d = dir(fullfile(projectRoot, 'data'));

% 1. Filter for items that are FOLDERS ([d.isdir])
% 2. Filter out self-references ('.', '..')
isActualFolder = [d.isdir] & ~ismember({d.name}, {'.', '..'});

% Extract the names of only the folders
folderNames = {d(isActualFolder).name};

% Final cleanup: remove .DS_Store (just in case it was flagged as a dir,
% though it is usually a file)
folderNames = folderNames(~strcmp(folderNames, '.DS_Store'));
folderNames = sort(folderNames);

end


function [cmd,cmd_fr] = get_cmd()
% GET_CMD  Generates shell commands to launch external MATLAB batches.
% This is used to solve SDP problems in a separate process to prevent
% "Out of Memory" or segmentation faults from crashing the main session.

projectRoot = get_project_root();
srcDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(projectRoot);
commonSrcDir = fullfile(packageRoot, 'src');
yalmipDir = fullfile(packageRoot, 'YALMIP-master');
projectRoot = strrep(projectRoot, '''', '''''');
srcDir = strrep(srcDir, '''', '''''');
commonSrcDir = strrep(commonSrcDir, '''', '''''');
yalmipDir = strrep(yalmipDir, '''', '''''');
matlabExecutable = fullfile(matlabroot, 'bin', 'matlab');

% Command for standard SDP relaxation
cmd = [
    matlabExecutable, ' -batch ', ...
    '"addpath(''', srcDir, '''); addpath(''', commonSrcDir, ...
    '''); addpath(genpath(''', yalmipDir, ''')); cd(''', projectRoot, ...
    '''); Robust_SAT_SDP" 2>/dev/null'
    ];

% Command for facial-reduced SDP relaxation
cmd_fr = [
    matlabExecutable, ' -batch ', ...
    '"addpath(''', srcDir, '''); addpath(''', commonSrcDir, ...
    '''); addpath(genpath(''', yalmipDir, ''')); cd(''', projectRoot, ...
    '''); Robust_SAT_SDP_FR" 2>/dev/null'
    ];


end


function projectRoot = get_project_root()
% GET_PROJECT_ROOT  Returns the SAT experiment directory.

srcDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(srcDir);

end



function output = save_instance(name,model,m,n)
% SAVE_INSTANCE  Packages basic problem information into a results structure.

output.name = regexprep(name, '\.cnf$', '');
output.n = n;
output.m = m;

if isfield(model,'obj')
    output.obj = model.obj;
end
if isfield(model,'Q')
    output.Q = model.Q;
end

end
