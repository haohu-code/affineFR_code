function BIPtable = BIP_test(maxsize,checkpointFile)
% this file demonstrates the algorithm described in
% Affine FR: an effective facial reduction algorithm for semidefinite
% relaxations of combinatorial problems
% https://arxiv.org/abs/2402.11796
% authors: Hao Hu and Boshi Yang
% last update: May-07-2024

if nargin < 1
    maxsize = 10000;
end
if nargin < 2
    checkpointFile = "";
end
checkpointFile = string(checkpointFile);

% Locate the MIPLIB data and the shared facial-reduction routines.
scriptDir = fileparts(mfilename('fullpath'));
miplibDir = fileparts(scriptDir);
packageRoot = fileparts(miplibDir);
dataDir = fullfile(miplibDir,'data');
addpath(fullfile(packageRoot,'src'));

% load the name of the mixed binary problems
% the problem is listed by the number of variables is increasing order
load(fullfile(miplibDir,'prob_list.mat'),'prob_list');

% MIPLIB files are third-party data and are not distributed with this code.
% Check the subset needed for this run before starting any computation.
requiredMask = cellfun(@(p) p.n <= maxsize,prob_list);
requiredNames = cellfun(@(p) p.name,prob_list(requiredMask), ...
    'UniformOutput',false);
installed = cellfun(@(name) isfile(fullfile(dataDir,name)),requiredNames);
if any(~installed)
    missingNames = requiredNames(~installed);
    preview = strjoin(missingNames(1:min(5,numel(missingNames))),', ');
    error('AffineFR:MIPLIBDataMissing', ...
        ['The MIPLIB instance files are not distributed with this code. ' ...
         '%d of %d files required for this run are missing from:\n  %s\n' ...
         'First missing files: %s\n' ...
         'See MIPLIB/README.md for download and installation instructions.'], ...
        nnz(~installed),numel(requiredNames),dataDir,preview);
end

% initialize the variables
np = length(prob_list); % number of problems
BIPtable = cell(np,1); % save the experiment results for table

% start the test
fprintf('BIP test at %s\n',string(datetime))
BIP_print([],false,true);
for i = 1:np
    % initialize placeholder for experimental results
    BIPout = [];
    
    % only solve instances with at most maxsize variables
    n = prob_list{i}.n;
    if n > maxsize
        break
    end

    % load the problem and make some corrections
    prob = gurobi_read(fullfile(dataDir,prob_list{i}.name));
    prob = prob_correction(prob,prob_list{i}.name(1:end-4));

    % save the problem information for latex printing
    BIPout.pinfo = prob_info(prob);

    % applying different FRA methods
    BIPout.FRA{1} = FRAcompute(prob,'affineFR');  % affine FR
    BIPout.FRA{2} = FRAcompute(prob,'partialFR_D');  % partial FR with D cone
    BIPout.FRA{3} = FRAcompute(prob,'partialFR_DD'); % partial FR with DD cone

    % print the result
    BIP_print(BIPout,false,false);

    % save information for latex table
    BIPtable{i} = BIP_statistics(BIPout);

    % Save a recoverable checkpoint after each completed instance.
    if strlength(checkpointFile) > 0
        save(checkpointFile,'BIPtable','maxsize','-v7.3');
    end
end

fprintf('The experiment is finished at %s\n',string(datetime))
fprintf('\n\n\n\n\n\n')


end

function BIPtable = BIP_statistics(BIPout)

% save instance information
BIPtable.pinfo = BIPout.pinfo;

% save Affine FR/Partial FR (D)/(DD)
for j = 1:3
    BIPtable.FRA{j}.fra_time = BIPout.FRA{j}.fra_time;
    BIPtable.FRA{j}.frv_time = BIPout.FRA{j}.frv_time;
    BIPtable.FRA{j}.total_time = BIPout.FRA{j}.total_time;
    BIPtable.FRA{j}.r = BIPout.FRA{j}.r;
    BIPtable.FRA{j}.status = BIPout.FRA{j}.status;
end

end

function BIP_print(results,is_tex,mode)
% INPUT:
% 1. results a structure containing the experiment results
% 2. is_tex do we need to add & sign for latex
% 3. mode = true if this is for printing header
%    mode = false if this is for printing experiment results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print the instance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mode
    myprint('%-23s','Instance',is_tex,mode)
    myprint('%-7s','Var',is_tex,mode)
else
    pinfo = results.pinfo;
    myprint('%-23s',pinfo.name,is_tex,~mode)
    myprint('%-7d',pinfo.n,is_tex,~mode)
end

myprint('%-1s','|',is_tex,~is_tex)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print FRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fratitle{1} = 'AFR';
fratitle{2} = 'PFR(D)';
fratitle{3} = 'PFR(DD)';
for j = 1:3
    if mode
        myprint('%-9s',fratitle{j},is_tex,mode)
        myprint('%-9s','ratio',is_tex,mode)
        myprint('%-6s','time',is_tex,mode)
    elseif isfield(results,'FRA')
        FRA = results.FRA{j};
        myprint('%-9d',FRA.r,is_tex,~mode)
        myprint('%-9.3f',FRA.r/(pinfo.n+1),is_tex,~mode)
        myprint('%-6.2f',FRA.total_time,is_tex,~mode) 
    else
        myprint('%-24s','Skip',is_tex,~mode)
    end
    myprint('%-1s','|',is_tex,~is_tex)
end

if is_tex
    fprintf('\\hline\n')
else
    fprintf('\n')
end

end

function myprint(s1,s2,istex,isprint)
% print s1 with argument s2 if flag = true

if isprint == true
    fprintf(s1,s2)
    if istex
        fprintf(' & ')
    end
end
end

function pinfo = prob_info(prob)
% save information about the problem instances for printing
pinfo.name = prob.modelname;        % model name
pinfo.n = length(prob.obj);         % number of variables
pinfo.nb = sum(or(prob.vtype=='I',prob.vtype=='B')); % binary variables
pinfo.m = length(prob.sense);       % total constraints
pinfo.m1 = sum(prob.sense=='=');    % equality constraints
pinfo.m2 = sum(prob.sense~='=');    % inequality constraints
pinfo.nf = sum(prob.lb == prob.ub); % number of fixed variables
pinfo.nl = sum(isfinite(prob.lb)); % number of lower bounds
pinfo.nu = sum(isfinite(prob.ub)); % number of lower bounds
pinfo.vtype = prob.vtype;    % variable type, continuous/binary/integer
end


function prob = prob_correction(prob,pname)

% add model name if necessary
if ~isfield(prob,'modelname')
    prob.modelname = pname;
end

% there are two bad models in the library
if strcmp(prob.modelname,'mas74')||strcmp(prob.modelname,'mas76')
    % the upper bound of variable x_{151} is 1e12
    % we scale the data for numerical stability
    prob.ub(151) = prob.ub(151)/1e6;
    prob.A(:,151) = prob.A(:,151)*1e6;
end

end
