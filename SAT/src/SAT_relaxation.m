% SAT_RELAXATION  Solves various convex relaxations (SDP, RLT) for the SAT ILP model.
%
%   output = SAT_relaxation(model, type, V) constructs and solves a relaxation
%   of the binary satisfiability problem using YALMIP and MOSEK.
%
%   Inputs:
%       model - Struct with fields A, rhs (constraints), obj (linear obj), Q (opt quadratic obj)
%       type  - String; 'SDP_SHOR', 'SDP_STRONG', or 'RLT'
%       V     - (Optional) Facial range vector for dimensionality reduction (Facial Reduction)
%
%   Outputs:
%       output - Struct containing objval, runtime, status, and YALMIP info.

function output = SAT_relaxation(model,type,V)

n = size(model.A,2);
% Lift the constraints: [A x >= b] becomes [Ahat * [1; x] >= 0] where Ahat = [-b, A]
Ahat = [-model.rhs model.A];

% Initialize symmetric matrix variable Y of size (n+1)x(n+1)
% Y is the lifted variable representing [1; x] * [1; x]^T
Y = sdpvar(n+1);

% --- Fundamental Constraints ---
% 1. Homogeneous: Y(1,1) == 1
% 2. Binary property (x_i^2 = x_i): implies diag(Y) == Y(:,1)
F = [Y(1,1) == 1; diag(Y) == Y(:,1)];

% --- Semidefinite Programming (SDP) Specific Preprocessing ---
if strcmp(type(1:3),'SDP')
    if nargin<3
        % Standard PSD constraint: Y must be positive semidefinite
        F = [F; Y>=0];
    else
        % Facial Reduction (FR): Y = V*R*V' where R is a smaller PSD matrix
        r = size(V,2);
        R = sdpvar(r);
        F = [F; Y == V*R*V'; R>=0];

        % Filter out redundant constraints relative to the range space of V
        idx = sum(abs(Ahat*V),2) < 1e-8;
        Ahat = Ahat(~idx,:);
    end
end

% --- Relaxation-Specific Constraints ---
if strcmp(type,'SDP_SHOR')
    % Classical Shor's SDP relaxation
    % Ensures that for each clause [1; x]^T * a_i >= 0, we have E[[1; x] * a_i^T]_diag >= 0
    F = [F; Ahat*diag(Y) >= 0];

elseif strcmp(type,'SDP_STRONG')
    % Lasserre-style / More robust relaxation
    % 1. Element-wise non-negativity (since x is binary)
    mask = triu(true(size(Y)), 1);
    F = [F; Y(mask) >= 0];

    % 2. First-order RLT constraints: Ahat * Y >= 0
    % This enforces Ahat * Y(:,j) >= 0 for all columns j, which is stronger than just Ahat * diag(Y) >= 0
    F = [F; Ahat*Y >=0];

elseif strcmp(type,'RLT')
    % Reformulation-Linearization Technique (Pure Linear approach)
    % Linear bounds: A * diag(Y) >= b
    F = [F; Ahat*diag(Y) >= 0];

    % McCormick Envelopes for variable products Y_ij = x_i * x_j
    mask = triu(true(size(Y)), 1);
    e = ones(n, 1);
    F = [F; Y(mask) >= 0; Y(mask) <= 1]; % 0 <= x_i*x_j <= 1

    % Linear inequalities defining the product of binary variables
    F = [F, get_lt(1 - Y(2:end,1)*e' - e*Y(2:end,1)' + Y(2:end, 2:end)) >= 0,...
        get_lt(Y(2:end,1)*e' - Y(2:end, 2:end)) >= 0];
else
    error('Invalid relaxation type: %s\n', type)
end



% --- Objective Function Selection ---
if isfield(model,'Q')
    % Quadratic objective: c^T*x + x^T*Q*x -> c^T*Y(2:end,1) + trace(Q*Y(2:end,2:end))
    obj = model.obj'*Y(2:end,1) + trace(blkdiag(0,model.Q)*Y);
else
    % Linear objective: c^T*x -> c^T*Y(2:end,1)
    obj = model.obj'*Y(2:end,1);
end

% Numerical settings for MOSEK
% opts = sdpsettings('solver','mosek','debug',1,'verbose',0);
opts = sdpsettings('solver','mosek','debug',1,'verbose',0,...
    'savesolveroutput',1); % Enable raw output capture
opts.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 3600; % Time limit in seconds


% Call YALMIP Optimizer
output = optimize(F,obj,opts);

% Organize and standardize output results
output.runtime = output.solvertime;
output.solver_message = output.info;
output.status = output.info;
output.objval = [];
if isfield(output,'solveroutput') && ...
        isstruct(output.solveroutput) && ...
        isfield(output.solveroutput,'res') && ...
        isstruct(output.solveroutput.res) && ...
        isfield(output.solveroutput.res,'sol') && ...
        isstruct(output.solveroutput.res.sol) && ...
        isfield(output.solveroutput.res.sol,'itr') && ...
        isstruct(output.solveroutput.res.sol.itr) && ...
        isfield(output.solveroutput.res.sol.itr,'dobjval')
    % Report Mosek's dual objective, consistent with the numerical tables.
    output.objval = -output.solveroutput.res.sol.itr.dobjval;
end

if and(strcmp(type(1:3),'SDP'),nargin==3)
    output.mr = sum(~idx); % Store count of reduced constraints
end


if contains(output.status,'Successfully solved')
    output.status = 'OPTIMAL';
elseif contains(output.status,'Infeasible problem')
    output.status = 'INFEASIBLE';
elseif contains(output.status,'Numerical problems')
    output.status = 'NUMERIC';
elseif contains(output.status, 'Unbounded objective function')
    output.status = 'UNBOUNDED';
elseif contains(output.status, 'Other identified error')
    output.status = 'ERROR';
else
    warning('AffineFR:UnknownSolverStatus', ...
        'Unrecognized solver status: %s',output.solver_message);
    output.status = 'ERROR';
end

end


function x = get_lt(X)
% return the lower triangular part of a matrix.
n = size(X);
idx = find(tril(ones(n)));
x = X(idx);

end
