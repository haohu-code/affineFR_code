% Run the 14 weighted-SAT experiments reported in manuscript Table 6.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
packageRoot = fileparts(projectRoot);
resultsDir = fullfile(projectRoot,'results');

addpath(scriptDir);
addpath(fullfile(projectRoot,'src'));
addpath(fullfile(packageRoot,'src'));
addpath(genpath(fullfile(packageRoot,'YALMIP-master')));
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

originalDirectory = pwd;
directoryCleanup = onCleanup(@() cd(originalDirectory));
cd(projectRoot);

instances = { ...
    'anomaly','blocksworld'; ...
    'medium','blocksworld'; ...
    'jnh14','jnh'; ...
    'jnh210','jnh'; ...
    'jnh215','jnh'; ...
    'jnh302','jnh'; ...
    'jnh309','jnh'; ...
    'jnh7','jnh'; ...
    'par8-1','parity'; ...
    'par8-2','parity'; ...
    'par8-3','parity'; ...
    'par8-4','parity'; ...
    'par8-5','parity'; ...
    'ssa0432-003','ssa'};

% These unreduced SDPs previously exhausted memory. They are deliberately
% not constructed or solved by this script.
unsafeOriginalSDP = { ...
    'par8-1'; ...
    'par8-2'; ...
    'par8-3'; ...
    'par8-4'; ...
    'par8-5'; ...
    'ssa0432-003'};

timestamp = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
rawResultsFile = fullfile(resultsDir, ...
    "table6_raw_" + timestamp + ".mat");

metadata = struct;
metadata.createdAt = char(datetime('now','TimeZone','local'));
metadata.matlabVersion = version;
metadata.matlabRelease = version('-release');
metadata.computer = computer;
metadata.gurobiFunction = which('gurobi');
metadata.mosekFunction = which('mosekopt');
metadata.objectiveFile = 'data/table6_objectives.mat';
metadata.relaxation = 'SDP_STRONG';
metadata.unsafeOriginalSDP = unsafeOriginalSDP;
metadata.unsafePolicy = [ ...
    'The listed unreduced SDPs are recorded as SKIPPED_OOM_SAFETY ', ...
    'and are not constructed or solved.'];

table6Results = cell(size(instances,1),1);
save(rawResultsFile,'table6Results','instances','metadata','-v7.3');

for instanceIndex = 1:size(instances,1)
    instanceName = instances{instanceIndex,1};
    datasetName = instances{instanceIndex,2};

    options = struct;
    options.mode = 2;
    options.obj_mode = 'linear';
    options.relaxation = 'SDP_STRONG';
    options.robust = true;
    options.computeV = true;
    options.datasets = {datasetName};
    options.instances = {instanceName};
    options.useTable6Objectives = true;
    options.skipOriginalSDPInstances = unsafeOriginalSDP;

    fprintf('\nTable 6 instance %d/%d: %s\n', ...
        instanceIndex,size(instances,1),instanceName);
    experimentOutput = SAT_batch_test(options);
    datasetOutput = experimentOutput{1};
    completed = find(~cellfun(@isempty,datasetOutput));
    if numel(completed) ~= 1
        error('Expected one result for %s, found %d.', ...
            instanceName,numel(completed));
    end

    table6Results{instanceIndex} = datasetOutput{completed};
    table6Results{instanceIndex}.runOptions = options;

    % Preserve every completed instance even if a later solver run fails.
    save(rawResultsFile,'table6Results','instances','metadata','-v7.3');
end

generatedFiles = generate_table6_outputs( ...
    table6Results,resultsDir,timestamp);

fprintf('\nGenerated Table 6 files:\n');
fprintf('  Raw results: %s\n',rawResultsFile);
fprintf('  CSV summary: %s\n',generatedFiles.csvFile);
fprintf('  MAT summary: %s\n',generatedFiles.matFile);
fprintf('  LaTeX table: %s\n',generatedFiles.latexFile);
