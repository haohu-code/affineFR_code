function compute_bejing_sparse_ranks
% Independent sparse-QR rank check for the six nonzero Bejing matrices.

root = fileparts(mfilename('fullpath'));
matrix_dir = fullfile(root,'results','bejing_rank_matrices');
files = dir(fullfile(matrix_dir,'*.mat'));
files = files(~[files.isdir]);

member = strings(length(files),1);
num_vars = zeros(length(files),1);
selected_rows = zeros(length(files),1);
structural_rank = zeros(length(files),1);
numerical_rank = zeros(length(files),1);
reduced_order = zeros(length(files),1);
qr_runtime_seconds = zeros(length(files),1);

for i = 1:length(files)
    data = load(fullfile(files(i).folder,files(i).name));
    N = data.N;
    rank_start = tic;
    [~,R,~] = qr(N,'vector');
    qr_runtime_seconds(i) = toc(rank_start);
    diagonal = full(abs(diag(R)));
    % FRAcompute treats eigenvalues of W=N'*N above 1e-5 as positive.
    % The corresponding singular-value threshold for N is sqrt(1e-5).
    numerical_rank(i) = sum(diagonal > sqrt(1e-5));
    member(i) = string(erase(files(i).name,'.mat')) + ".cnf";
    num_vars(i) = data.n;
    selected_rows(i) = size(N,1);
    structural_rank(i) = sprank(N);
    reduced_order(i) = data.n + 1 - numerical_rank(i);
    fprintf('%s: rank=%d, r=%d, QR time=%.3f s\n', ...
        member(i),numerical_rank(i),reduced_order(i),qr_runtime_seconds(i));
end

results = table(member,num_vars,selected_rows,structural_rank, ...
    numerical_rank,reduced_order,qr_runtime_seconds);
writetable(results,fullfile(root,'results','bejing_sparse_ranks.csv'));
end
