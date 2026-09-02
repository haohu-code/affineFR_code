function PFRout = partialFR(prob,cone_type)
% implementation of partial FR with D cone and DD cone
fra_time = tic;

% load the parameters
if isfield(prob,'obj')
    n = length(prob.obj);
elseif isfield(prob,'vtype')
    n = length(prob.vtype);
end

% the indices of the lower triangular entries of the matrix variable
ridx = offDiagidx(n+1);

% get the matrix representation for the dual problem
%m1 and m2 is number of dual variables associated with equality and
%inequality constraints
[Aadj,m1,m2] = get_Adj(prob,n);
Aadj = Aadj(ridx,:);
m = m1+m2;

% get the cone
D = get_cone(n+1,cone_type);
D = D(ridx,:);

% apply partial FR
output = Pfr_opt(Aadj,D,m1,m2,n);

% obtain the exposing vector
if strcmp(output.status,'OPTIMAL')
    % get the exposing vector
    y = output.x(1:m); % dual variable except arrow
    z = output.x(m+1:m+n); % dual variable associated with arrow
    W = offDiag2M(Aadj*[y;z],n+1); % exposing vector
    W = (W+W')/2;

    d = output.x(m+n+1:m+n+size(D,2)); % cone variables
    t = output.x(m+n+size(D,2)+1:end); % auxiliary variables
    PFRout.ymax = max(d+t);

else
    % if the problem is not solved correctly
    W = sparse(n+1,n+1);
    PFRout.r = n + 1;
    PFRout.ymax = 0;
end

% save the results;
PFRout.W = W;
PFRout.solver_output = output;
PFRout.status = output.status;
PFRout.fra_time = toc(fra_time); % total time

end


function output = Pfr_opt(Aadj,D,m1,m2,n)
% apply partial FR
% Aadj the matrix representation of the dual constraint
% D the matrix presentation of the cone (D) or (DD)
% m1 number of dual var for equalities
% m2 number of dual var for inequalities
% n number of dual for arrow constraints

m = m1 + m2;
[k,nd] = size(D);
nv = m+n+2*nd;
prob.A = sparse([Aadj -D sparse(k,nd); ...
    sparse(nd,m+n) speye(nd) -speye(nd)]);
k = size(prob.A,1);
prob.rhs = zeros(k,1);
prob.sense = [char(ones(size(Aadj,1),1)*'='); char(ones(nd,1)*'>')];
prob.ub = Inf*ones(nv,1);
prob.ub(end-nd+1:end) = 1;
prob.lb = zeros(nv,1);
prob.lb(1:m1) = -Inf; % the dual var for equalities are free
prob.lb(m+1:m+n) = -Inf; % the dual var for arrow constraints are free
prob.vtype = char(ones(nv,1)*'C');
prob.obj = zeros(nv,1);
prob.obj(end-nd+1:end) = 1;
prob.modelsense = 'max';

% keyboard
pars = [];
pars.outputflag = 0; % do not print
pars.TimeLimit = 60*5; % we should never spend too much time in preprocessing
output = gurobi(prob,pars);

% If the problem is infeasible, then we try to solve the problem again with
% numericfocus = 3. This ensures that the infeasibility is not due to
% bad models.
if strcmp(output.status,'INFEASIBLE')
    pars.NumericFocus = 3;
    output = gurobi(prob,pars);
end

% save the solution
if strcmp(output.status,'OPTIMAL')
    x = output.x;
    output.y = x(1:m);
    output.z = x(m+1:m+n);
    output.d = x(m+n+1:end);
    output.prob = prob;
    output.pres = norm(prob.A*x - prob.rhs,2) + sum(max(-x(m1+1:m),0))+ sum(max(-x(m+n+1:end),0));
else
    output.pres = -1;
end

end

function [Aadj,m1,m2] = get_Adj(prob,n)
% get the linear constraint for the dual problem
[Aineq,bineq,Aeq,beq] = FRAformat(prob);
m1 = length(beq); % associated with equality
m2 = length(bineq); % associated with inequality

% the matrix representation of the constraints
At = [[-beq Aeq]' [-bineq Aineq]'];

% the matrix representation of the shor's data matrices
MS = shor_mat(n);

% the matrix representation of the arrow constraints
MA = arrow_adj(prob);

% combine everything A^{*} (A adjoint)
Aadj = [MS*At -MA];
end

function M = shor_mat(n)
% M:Rn->Sn+1 is the matrix representation of linear constraint in the
% Shor's SDP relaxation
M = zeros(n+1);
M(1:end,1) = 1;
idx1 = find(M);

M = zeros(n+1);
M(1,2:end) = 1;
idx2 = find(M);

I = [idx1;idx2];
J = [1:n+1 2:n+1]';
K = [1; ones(2*n,1)/2];
M = sparse(I,J,K,(n+1)^2,n+1);
end

function M = arrow_adj(prob)
% For BIP with n variables and the set of binary variables Bidx
% M:Rn->Sn+1 is the matrix representation of the adjoint of the arrrow operator
% For FR for QCQP, the first entry is zero and thus assumed to be removed here.

% obtain indices
n = length(prob.vtype);
Bidx = find(or(prob.vtype=='I',prob.vtype=='B'))+1;
BidxM = sub2ind([n+1 n+1],Bidx,Bidx);
nb = length(Bidx); % the number of binary variables

% the first column
M = zeros(n+1);
M(Bidx,1) = 1;
idx1 = find(M);

% the first row
M = zeros(n+1);
M(1,Bidx) = 1;
idx2 = find(M);

% the diagonal entries
M = zeros(n+1);
M(BidxM) = 1;
idx3 = find(M);

I = [idx1;idx2;idx3];
J = kron(ones(3,1),(1:nb)');
K = [-ones(nb,1);-ones(nb,1);2*ones(nb,1)];
M = sparse(I,J,K,(n+1)^2,n);

end

function D = get_cone(n,cone_type)
% get the extreme rays of the cone
switch cone_type
    case 'D'
        % non-negative diagonal cone
        D = getDcone(n);
    case 'DD'
        % get the extreme rays of diagonally dominant cone
        % D = getDDcone(n);
        D = getDDcone_eco(n);
end

end

function D = getDcone(n)
% Return the set of extreme rays of the non-negative diagonally cone of order n.
% Each extreme ray is represented as a column vector in D.

% get the indices of pairs
I = find(eye(n));
J = [1:n]';
K = ones(n,1);
D = sparse(I,J,K,n^2,n);
end


function D = getDDcone(n)
% Return the set of extreme rays of the diagonally dominant cone of order n.
% Each extreme ray is represented as a column vector in D.
% this version is a direct implementation and is inefficient,
% we leave it here for test and verification purpose.
% for an efficient implementation, please check ``D = getDDcone_eco(n)''

m = n*(n-1); % the number of off-diagonals

% get the indices of paris
I = kron([1:n]',ones(n,1));
J = kron(ones(n,1),[1:n]');
IJ = [min([I J],[],2) max([I J],[],2)];
IJ = unique(IJ,'rows');
IJ(IJ(:,1)==IJ(:,2),:) = [];
I = IJ(:,1);
J = IJ(:,2);
I = sub2ind([n n],[I;J;I;J],[I;I;J;J]);
I = reshape(I,m/2,4)';
I = [[1:n+1:n^2]';I(:);I(:)]';
J = [[1:n]'; kron([1:m/2]'+n,ones(4,1));kron([1:m/2]'+n+m/2,ones(4,1))];
K = [ones(n,1);ones(4*m/2,1); kron(ones(m/2,1),[1 -1 -1 1]')];
D = sparse(I,J,K,n^2,n^2);
end


function D = getDDcone_eco(n)
% Return the set of extreme rays of the diagonally dominant cone of order n.
% Each extreme ray is represented as a column vector in D.
%
% eco means:
% the set of extreme rays is quadratic in the number of variables
% it becomes intractable very soon as n becomes a few thousands
% we use the arrow structure to remove some redundant extreme rays here

% get the indices of pairs
I = ones(n-1,1);
J = [2:n]';
I = sub2ind([n n],[I;J;I;J],[I;I;J;J]);
I = reshape(I,n-1,4)';

I = [[1:n+1:n^2]';I(:);I(:)];
J = [[1:n]'; kron([1:n-1]'+n,ones(4,1));kron([1:n-1]'+2*n-1,ones(4,1))];
K = [ones(n,1);ones(4*(n-1),1); kron(ones(n-1,1),[1 -1 -1 1]')];
D = sparse(I,J,K,n^2,n+2*(n-1));
end


function idx = offDiagidx(n)
% get the linear indices of the off-diagonal entries of n by n matrix
M = zeros(n);
M(1:end) = 1:n^2;
M = tril(M);
idx = M(tril(true(n)));
end

function M = offDiag2M(z,n)
% given a vector z associated with the lower triangular part of the matrix
% M, this function recovers M
M = zeros(n);
idx = offDiagidx(n);
M(idx) = z;
M = M + tril(M,-1)';
end
