function model = SAT_ILP_FEA(satProblem,obj_mode)
% This function returns the ILP formulation for the SAT problem
% Constraints: Building the matrix A and RHS b
% Rearranging the formula: sum(y_positive) + sum(1 - y_negative) >= z_j
% Becomes: sum(y_positive) - sum(y_negative) - z_j >= - (count of negative literals)


n = satProblem.numVars;      % Number of y variables
m = satProblem.numClauses;   % Number of z variables
m = length(satProblem.clauses);

% --- Gurobi Model Construction ---

% 1. Variable Types: 'B' for Binary {0,1}
% First n variables are y_i, next m variables are z_j
model.vtype = repmat('B', 1, n);
model.lb = zeros(n,1);
model.ub = ones(n,1);
model.obj = zeros(n,1);
% assign a random weight
if strcmp(obj_mode,'linear')
    model.obj = round(abs(randn(n,1))*10);
elseif strcmp(obj_mode,'quadratic')
    % Q = abs(randn(n));
    Q = round(abs(randn(n))*10);
    Q = tril(Q)+tril(Q,-1)';
    model.Q = sparse(Q);
end
model.modelsense = 'min';

% 2. Preallocate for Speed
% Count total literals across all clauses to get exact size
numEntries = sum(cellfun(@length, satProblem.clauses));

rowIdx = zeros(numEntries, 1);
colIdx = zeros(numEntries, 1);
valIdx = zeros(numEntries, 1);
counter = 0; % Points to the current position in our preallocated vectors
b = zeros(m, 1);
for i = 1:m
    currentClause = satProblem.clauses{i}; % Access the cell
    for j = 1:length(currentClause)
        lit = currentClause(j);
        counter = counter + 1;
        rowIdx(counter) = i;
        colIdx(counter) = abs(lit);
        valIdx(counter) = sign(lit); % Shortcut: 1 if positive, -1 if negative
    end
    % The constant term from (1 - y_i) moves to the RHS
    b(i) = 1 - sum(currentClause < 0);
end


model.A = sparse(rowIdx, colIdx, valIdx, m, n);
model.rhs = b;
model.sense = repmat('>', m, 1); % All constraints are >=


end