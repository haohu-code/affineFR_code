function [Aineq,bineq,Aeq,beq] = FRAformat(prob)
% convert constraints in gurobi format into the form
% Aineq*x <= bineq and Aeq*x = beq

% convert the instance
eq_idx = prob.sense=='=';
leq_idx = prob.sense=='<';
geq_idx = prob.sense=='>';
Aeq = prob.A(eq_idx,:);
beq = prob.rhs(eq_idx,:);
Aineq = [prob.A(leq_idx,:); -prob.A(geq_idx,:)];
bineq = [prob.rhs(leq_idx,:); -prob.rhs(geq_idx,:)];

% lower and upper bounds as inequality
n = length(prob.obj);
lidx = ~isinf(prob.lb);
uidx = ~isinf(prob.ub);
L = speye(n);
L = L(lidx,:);
U = speye(n);
U = U(uidx,:);
l = prob.lb(lidx);
u = prob.ub(uidx);

Aineq = [Aineq;-L;U];
bineq = [bineq;-l;u];

end