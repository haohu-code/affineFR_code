function FRAout = FRAcompute(prob,method,computeV)
% compute an exposing vector and facial range vector for any SDP relaxation
% of the mix binary programming problem given in gurobi format
% INPUT: prob     - input problem in gurobi format
%        method   - facial reduction method
%        computeV - whether to construct a facial range vector (default true)
% OUTPUT: FRAout is a structure containing the following fields
%         Part 1:
%         W - the exposing vector, or [] when it is not constructed
%         solver_output - the solver output
%         status - the solver status
%         fra_time - the computational time of the FRA method
%         Part 2:
%         r - the order of the matrix variable after reduction
%         V - the facial range vector
%         frv_time - the time for computing V (part of SDP reformulation)
%         rank_time - the sparse-rank time when computeV is false

if nargin < 3
    computeV = true;
end

% Part 1: applying a special FRA to obtain an exposing vector
if strcmp(method,'affineFR')
    FRAout = affineFR(prob,computeV);  % affine FR
elseif strcmp(method,'partialFR_D')
    FRAout = partialFR(prob,'D');  % partial FR with D cone
elseif strcmp(method,'partialFR_DD')
    FRAout = partialFR(prob,'DD'); % partial FR with DD cone
end

if computeV
    % Part 2: construct a facial range vector for SDP reformulation.
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
    FRAout.rank_time = 0;
else
    % Rank-only mode for large instances. For affine FR, W=N'*N, so the
    % eigenvalue threshold 1e-5 corresponds to the singular-value threshold
    % sqrt(1e-5) for N. Sparse QR avoids forming W and does not construct V.
    rank_time = tic;
    if strcmp(method,'affineFR')
        [~,R,~] = qr(FRAout.N,'vector');
        rankW = sum(abs(diag(R)) > sqrt(1e-5));
        FRAout.r = size(FRAout.N,2)-rankW;
    else
        [~,R,~] = qr(FRAout.W,'vector');
        rankW = sum(abs(diag(R)) > 1e-5);
        FRAout.r = size(FRAout.W,1)-rankW;
    end
    FRAout.V = [];
    FRAout.frv_time = 0;
    FRAout.rank_time = toc(rank_time);
end

FRAout.computeV = logical(computeV);
FRAout.total_time = FRAout.fra_time + FRAout.frv_time + FRAout.rank_time;
end
