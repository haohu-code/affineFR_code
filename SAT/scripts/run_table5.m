% Reproduce Table 5.

% Locate the SAT-specific and shared source directories.
scriptDir = fileparts(mfilename('fullpath'));
satCodeDir = fileparts(scriptDir);
packageRoot = fileparts(satCodeDir);
addpath(fullfile(satCodeDir,'src'));
addpath(fullfile(packageRoot,'src'));

% Run affine FR on exactly the ten families reported in Table 5.
% Computing V is unnecessary because the table reports only reduced orders.
rng("default")
options.mode = 1;
options.robust = false;
options.computeV = false;
options.datasets = table5_families();

% Compute the per-instance results and the family-level averages.
[EXPoutputs,folderNames] = SAT_batch_test(options);
summary = summarize_table5(EXPoutputs,folderNames);
disp(summary)

% Save both the complete results and the summary used in the paper.
resultsDir = fullfile(satCodeDir,'results');
if ~isfolder(resultsDir)
    mkdir(resultsDir)
end
timestamp = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
save(fullfile(resultsDir,"table5_" + timestamp + ".mat"), ...
    'EXPoutputs','folderNames','options','summary','-v7.3');
writetable(summary,fullfile(resultsDir,"table5_" + timestamp + ".csv"));
