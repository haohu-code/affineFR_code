function AFRout = affineFR(prob,computeW)
% This function is an implementation of affine FR
% Find the affine hull of the polyhedron {x: Aineq*x<=bineq, Aeq*x=beq}
% PRIMAL max{0|Aineq*x<=bineq, Aeq*x=beq}
% DUAL   min{bineq'y+beq'z|Aineq'y+Aeq'z=0,y>=0}
% we will identify a dual optimal solution with the maximum number of
% nonzeros
%
% Affine FR: an effective facial reduction algorithm for semidefinite
% relaxations of combinatorial problems
% https://arxiv.org/abs/2402.11796
% authors: Hao Hu and Boshi Yang
% last update: May-07-2024

if nargin < 2
    computeW = true;
end

% setting
fra_time = tic;

% load the parameters
n = length(prob.vtype);

% get the data matrices
[Aineq,bineq,Aeq,beq] = FRAformat(prob);
Aeq1 = [Aineq'; bineq'];
Aeq2 = [Aeq'; beq'];

% find the affine hull
solvername = 'gurobi';
[y,output] = split_gurobi(Aeq1,Aeq2);


% if the affine hull is computed correctly
if output.flag
    idx = y(1:size(Aineq,1)) > 1/2;
    
    B = [Aeq;Aineq(idx,:)];
    d = [beq;bineq(idx,:)];
    AFRout.N = [-d B];
    if computeW
        W = AFRout.N'*AFRout.N;
        AFRout.W = (W+W')/2;
    else
        AFRout.W = [];
    end
    AFRout.flag = true;
    AFRout.ymax = max(y(1:size(Aineq,1)));
else
    % if the solver can't compute the affine hull
    AFRout.N = sparse(0,n+1);
    if computeW
        AFRout.W = sparse(n+1,n+1);
    else
        AFRout.W = [];
    end
    AFRout.flag = false;
    AFRout.ymax = 0;
end

% save the results
AFRout.solvername = solvername;
AFRout.solver_output = output;
AFRout.status = output.status;
AFRout.fra_time = toc(fra_time); % total time

end

function [y,output] = split_gurobi(A,B)
% find an element in the polyhedron
% P = {x | A*x + B*z = 0, x>=0,z free}
% such that x has the maximum number of non-zeros
%
% this is achieved by solving
% P = max{1^{T}u | A*u + Av + B*z = 0, u,v>=0, u<=1, z free}
% then x = u + v has the maximum number of non-zeros

% find implicit equalities
[m,n] = size(A);
[~,n2] = size(B);

% define the problem
prob.A = sparse([A A B]);
prob.rhs = zeros(m,1);
prob.sense = char(ones(m,1)*'=');
prob.lb = [zeros(2*n,1);-Inf*ones(n2,1)];
prob.ub = [ones(n,1);Inf*ones(n+n2,1)];
prob.obj = [ones(n,1);zeros(n+n2,1)];
prob.modelsense = 'max';

% solve the problem
pars = [];
pars.outputflag = 0;
pars.TimeLimit = 60*5; % we should never spend too much time in preprocessing
output = gurobi(prob,pars);

% If the problem is infeasible, then we try to solve the problem again with
% numericfocus = 3. This ensures that the infeasibility is not due to
% bad models.
if strcmp(output.status,'INFEASIBLE')
    pars.NumericFocus = 3;
    output = gurobi(prob,pars);
end

if (strcmp(output.status,'OPTIMAL')||strcmp(output.status,'SUBOPTIMAL'))
    % solved normally
    y = output.x;

    % Interior-point solution.
    % compute residual
    output.pres = sqrt(norm(prob.A*y,2)^2 + ...
        norm(min(y(1:2*n),0),2)^2 + ...
        norm(max(y(1:n)-1,0),2)^2);
    y = [y(1:n)+y(n+1:2*n); y(2*n+1:end)];
    output.y = y;
    % if residual is too big, then we do not consider it as reliable
    if output.pres < 1e-5
        output.flag = true;
    else
        output.flag = false;
        output.status = 'bigres';
    end
else
    y = [];
    output.pres = -1;
    output.flag = false;
end

end
