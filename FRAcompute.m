function FRAout = FRAcompute(prob,method)
% compute an exposing vector and facial range vector for any SDP relaxation
% of the mix binary programming problem given in gurobi format
% INPUT: prob - input problem in gurobi format
% OUTPUT: FRAout is a structure containing the following fields
%         Part 1:
%         W - the exposing vector obtained from the specified method
%         solver_output - the solver output
%         status - the solver status
%         fra_time - the computational time of the FRA method
%         Part 2:
%         r - the order of the matrix variable after reduction
%         V - the facial range vector
%         frv_time - the time for computing V (part of SDP reformulation)

% Part 1: applying a special FRA to obtain an exposing vector
if strcmp(method,'affineFR')
    FRAout = affineFR(prob);  % affine FR
elseif strcmp(method,'partialFR_D')
    FRAout = partialFR(prob,'D');  % partial FR with D cone
elseif strcmp(method,'partialFR_DD')
    FRAout = partialFR(prob,'DD'); % partial FR with DD cone
end

% Part 2: compute facial range vector for reformulating the SDP later
frv_time = tic;
try
    [V,d] = eig(full(FRAout.W),'vector');
catch ME
    warning(ME.message);
    fprintf('Warning: No reduction as the eigenvalue computation failed\n')
    n = length(FRAout.W);
    V = speye(n,n);
    d = zeros(n,1);
end

idx = d>1e-5;
FRAout.r = length(d)-sum(idx);
FRAout.V = V(:,~idx);
FRAout.frv_time = toc(frv_time);
FRAout.total_time = FRAout.frv_time + FRAout.fra_time;
end