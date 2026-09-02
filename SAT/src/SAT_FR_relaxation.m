% SAT_FR_RELAXATION  Solves the facially reduced convex relaxation for the SAT ILP model.
%
%   output = SAT_FR_relaxation(model, type, N) constructs and solves a facially
%   reduced SDP relaxation of the binary satisfiability problem using YALMIP and MOSEK.
%
%   Inputs:
%       model - Struct with fields A, rhs (constraints), obj (linear obj), Q (opt quadratic obj)
%       type  - String; (unused in current implementation, kept for signature compatibility)
%       N     - Normal subspace matrix determined from Facial Reduction Analysis (FRA)
%
%   Outputs:
%       output - Struct containing objval, runtime, status, vtime, mr, and YALMIP info.

function output = SAT_FR_relaxation(model,type,N)

% compute a sparse facial range vector
vtime = tic;
V = null(full(N),'r'); 
vtime = toc(vtime);


n = size(model.A,2);
% Lift the constraints: [A x >= b] becomes [Ahat * [1; x] >= 0] where Ahat = [-b, A]
Ahat = [-model.rhs model.A];


% Filter out redundant constraints relative to the range space of V
idx = sum(abs(Ahat*V),2) < 1e-8;
Ahat = Ahat(~idx,:);

% SDP variable
r = size(V,2);
R = sdpvar(r);
% add constraints
% psd and arrow
F = [R>=0;V(1,:)*R*V(1,:)' == 1; diag(V*R*V') == V*R*V(1,:)'];
% RLT constraints
F = [F; (Ahat*V)*R*V' >=0];
% non-negativity
F = [F; reshape(V*R*V',(n+1)^2,1)>=0];

% --- Objective Function Selection ---
if isfield(model,'Q')
    % Quadratic objective: c^T*x + x^T*Q*x -> c^T*Y(2:end,1) + trace(Q*Y(2:end,2:end))
    obj =  (model.obj'*V(2:end,:))*R*V(1,:)' + trace((V'*blkdiag(0,model.Q)*V)*R);
else
    % Linear objective: c^T*x -> c^T*Y(2:end,1)
    obj = (model.obj'*V(2:end,:))*R*V(1,:)';
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

output.vtime = vtime;
output.mr = sum(~idx); % Store count of reduced constraints



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
