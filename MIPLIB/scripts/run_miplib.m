% Run the complete MIPLIB preprocessing experiment and generate tables.

scriptDir = fileparts(mfilename('fullpath'));
miplibDir = fileparts(scriptDir);
packageRoot = fileparts(miplibDir);
resultsDir = fullfile(miplibDir,'results');

addpath(scriptDir);
addpath(fullfile(packageRoot,'src'));

if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

timestamp = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
rawResultsFile = fullfile(resultsDir,"miplib_raw_" + timestamp + ".mat");
maxsize = 10000;

metadata = struct;
metadata.createdAt = char(datetime('now','TimeZone','local'));
metadata.matlabVersion = version;
metadata.matlabRelease = version('-release');
metadata.computer = computer;
metadata.gurobiFunction = which('gurobi');
metadata.maxVariables = maxsize;
metadata.lpTimeLimitSeconds = 300;
metadata.eigenvalueThreshold = 1e-5;
metadata.constructFacialRangeVector = true;

% BIP_test validates the separately installed MIPLIB data before starting and
% writes a recoverable checkpoint after every completed instance.
BIPtable = BIP_test(maxsize,rawResultsFile);
save(rawResultsFile,'BIPtable','metadata','maxsize','-v7.3');

generatedFiles = generate_miplib_tables(BIPtable,resultsDir,timestamp);

fprintf('\nGenerated MIPLIB files:\n');
fprintf('  Raw results:       %s\n',rawResultsFile);
fprintf('  Reduction summary: %s\n',generatedFiles.reductionCSV);
fprintf('  Timing summary:    %s\n',generatedFiles.timingCSV);
fprintf('  Instance results:  %s\n',generatedFiles.instanceCSV);
fprintf('  Summary MAT file:  %s\n',generatedFiles.summaryMAT);
fprintf('  LaTeX tables:      %s\n',generatedFiles.latexFile);
