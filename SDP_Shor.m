function output = SDP_Shor(prob)
% this functions solve the Shor's SDP relaxation for a given BIP instance
% INPUT: prob is an BIP instance in gurobi format
% OUTPUT: a lower bound from the Shor's SDP relaxation

% load the setting
yalmip('clear')

% convert the problem into Aineq*x<=bineq and Aeq*x = beq format
n = length(prob.obj);
bidx = find(or(prob.vtype=='I',prob.vtype=='B'));
bidxM = sub2ind([n n],bidx,bidx); % the position of binary variables in the diagonal of the matrix variable X

% get the linear constraint for the dual problem
[Aineq,bineq,Aeq,beq] = FRAformat(prob);

% construct yalmip model
X = sdpvar(n);
x = sdpvar(n,1);
F = [Aeq*x == beq; Aineq*x<=bineq; X(bidxM) == x(bidx); [1 x';x X]>=0];
obj = prob.obj'*x;
opts = sdpsettings('solver','mosek','debug',1,'verbose',0,'cachesolvers',1);
opts.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 3600;
output = optimize(F,obj,opts);

x = value(x);
X = value(X);
lb = value(obj);
output.Y = [1 x';x X];
output.lb = lb;
end