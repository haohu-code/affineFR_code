function output = SDP_Shor_FR(prob,V)
% this functions solve the facially reduced Shor's SDP relaxation
% for a given BIP instance
% INPUT: prob is an BIP instance in gurobi format
%        V - the facial range vector
% OUTPUT: a lower bound from the Shor's SDP relaxation

% load the setting
yalmip('clear')
% infeasibility is detected without the need of solving SDP
if size(V,2) == 0
    output.lb = -Inf;
    output.time = 0;
    return
end

n = length(prob.vtype);
mode = 2;
% mode = 1 naive implementation for verification
% mode = 2 official implementation
if mode == 1
    % the index of the binary variables in terms of n by n matrix variable X
    bidx = find(or(prob.vtype=='I',prob.vtype=='B'));
    bidxM = sub2ind([n n],bidx,bidx);
    % get the linear constraint for the dual problem
    [Aineq,bineq,Aeq,beq] = FRAformat(prob);

    c = prob.obj;

    r = size(V,2);
    X = sdpvar(n);
    R = sdpvar(r);
    x = sdpvar(n,1);
    F = [Aeq*x == beq;
        Aineq*x<=bineq;
        X(bidxM) == x(bidx);
        [1 x';x X]  == V*R*V';
        R>=0];
    obj = c'*x;

else
    % the index of the binary variables in terms of n+1 matrix variable Y
    bidx = find(or(prob.vtype=='I',prob.vtype=='B'))+1;
    bidxM = sub2ind([n+1 n+1],bidx,bidx);

    % get the linear constraint for the dual problem
    [Aineq,bineq,Aeq,beq] = FRAformat(prob);

    c = prob.obj;

    r = size(V,2);
    R = sdpvar(r);
    F = [V(1,:)*R*V(1,:)' == 1;
        ([-beq Aeq]*V)*R*V(1,:)' == 0;
        ([-bineq Aineq]*V)*R*V(1,:)'<=0;
        diag_idx(V*R*V',bidxM) == V(bidx,:)*R*V(1,:)';
        R>=0];
    obj = c'*V(2:end,:)*R*V(1,:)';
end

% solve the problem
opts = sdpsettings('solver','mosek','debug',1,'verbose',0,'cachesolvers',1);
% opts.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 3600;
output = optimize(F,obj,opts);

% save the output
R = value(R);
lb = value(obj);
output.V = V;
output.R = R;
output.lb = lb;
output.time = output.solvertime;

end

function M2 = diag_idx(M,idx)
M2 = M(idx);
end