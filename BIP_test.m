function BIPtable = BIP_test(maxsize,solveSDP)
% this file demonstrates the algorithm described in
% Affine FR: an effective facial reduction algorithm for semidefinite
% relaxations of combinatorial problems
% https://arxiv.org/abs/2402.11796
% authors: Hao Hu and Boshi Yang
% last update: May-07-2024

if nargin < 1
    maxsize = 10000;
    solveSDP = false;
end

% add the path to the solvers in the linux server for MPC
addpath(genpath('/opt/gurobi911'))
addpath(genpath('/opt/mosek'))

% Add the current folder and its all subfolders to the path.
file_path = mfilename('fullpath');
idx = find(file_path == '/');
folder = file_path(1:idx(end)); % location of current m file
addpath(genpath(folder)); % add to path

% load the name of the mixed binary problems
% the problem is listed by the number of variables is increasing order
load('prob_list.mat','prob_list');

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
    try
        prob = gurobi_read([folder 'data/' prob_list{i}.name]);
    catch
        fprintf('this problem is not in data folder\n')
        continue
    end
    prob = prob_correction(prob,prob_list{i}.name(1:end-4));

    % save the problem information for latex printing
    BIPout.pinfo = prob_info(prob);

    % applying different FRA methods
    BIPout.FRA{1} = FRAcompute(prob,'affineFR');  % affine FR
    BIPout.FRA{2} = FRAcompute(prob,'partialFR_D');  % partial FR with D cone
    BIPout.FRA{3} = FRAcompute(prob,'partialFR_DD'); % partial FR with DD cone

    % solve (facially reduced) Shor's relaxation
    if solveSDP == true
        flag = false;
        for j = 1:3
            BIPout.Shor{j}.flag = false;
            % if there exists a postive reduction
            if BIPout.FRA{j}.r < n + 1
                % solve (facially reduced) Shor's relaxation
                BIPout.Shor{j} = SDP_Shor_FR(prob,BIPout.FRA{j}.V);
                BIPout.Shor{j}.flag = true;
                flag = true;
            end
        end
        % if there exists a positive reduction for any of the 3 methods
        BIPout.Shor{4}.flag = false;
        if flag == true
            % solve the Shor's SDP relaxation
            BIPout.Shor{4} = SDP_Shor(prob);
            BIPout.Shor{4}.flag = true;
        end
    end

    % print the result
    BIP_print(BIPout,false,false);

    % save information for latex table
    BIPtable{i} = BIP_statistics(BIPout,solveSDP);
end

fprintf('The experiment is finished at %s\n',string(datetime))
fprintf('\n\n\n\n\n\n')


end

function BIPtable = BIP_statistics(BIPout,solveSDP)

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

% save SDP results, if we solved any
if solveSDP == true
    for j = 1:4
        % save the flag if we solved the SDP
        BIPtable.Shor{j}.flag = BIPout.Shor{j}.flag;
        % if so, saved the details
        if BIPout.Shor{j}.flag
            BIPtable.Shor{j}.yalmiptime = BIPout.Shor{j}.yalmiptime;
            BIPtable.Shor{j}.solvertime = BIPout.Shor{j}.solvertime;
            BIPtable.Shor{j}.lb = BIPout.Shor{j}.lb;
        elseif BIPout.Shor{4}.flag
            % if one of the FR has no reduction, then save Shor's
            % result
            BIPtable.Shor{j}.yalmiptime = BIPout.Shor{4}.yalmiptime;
            BIPtable.Shor{j}.solvertime = BIPout.Shor{4}.solvertime;
            BIPtable.Shor{j}.lb = BIPout.Shor{4}.lb;
        end
    end
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print SDP Shor's relaxation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shortitle{1} = 'Shor_AFR_lb';
shortitle{2} = 'Shor_PFR(D)_lb';
shortitle{3} = 'Shor_PFR(DD)_lb';
shortitle{4} = 'Shor_lb';
for j = 1:4
    if mode
        myprint('%-16s',shortitle{j},is_tex,mode)
        myprint('%-8s','time',is_tex,mode)
    elseif isfield(results,'Shor')&&results.Shor{j}.flag
        Shor = results.Shor{j};
        myprint('%-16.2e',Shor.lb,is_tex,~mode)
        myprint('%-8.2f',Shor.solvertime,is_tex,~mode)
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