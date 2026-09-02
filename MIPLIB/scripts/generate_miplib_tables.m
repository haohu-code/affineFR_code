function files = generate_miplib_tables(BIPtable,resultsDir,timestamp)
%GENERATE_MIPLIB_TABLES Aggregate MIPLIB results and write CSV/LaTeX files.

if nargin < 2 || strlength(string(resultsDir)) == 0
    scriptDir = fileparts(mfilename('fullpath'));
    resultsDir = fullfile(fileparts(scriptDir),'results');
end
if nargin < 3 || strlength(string(timestamp)) == 0
    timestamp = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
end
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

completed = ~cellfun(@isempty,BIPtable);
BIPtable = BIPtable(completed);
if isempty(BIPtable)
    error('No completed MIPLIB instances were supplied.');
end

methodNames = ["Affine FR"; "Partial FR (D)"; "Partial FR (DD)"];
latexMethodNames = {'Affine FR', ...
    'Partial FR ($\mathcal{D}$)', ...
    'Partial FR ($\mathcal{DD}$)'};
numMethods = numel(methodNames);
numInstances = numel(BIPtable);

instanceNames = strings(numInstances,1);
numVariables = zeros(numInstances,1);
originalOrder = zeros(numInstances,1);
reducedOrder = zeros(numInstances,numMethods);
orderRatio = zeros(numInstances,numMethods);
auxiliaryLPTime = zeros(numInstances,numMethods);
basisTime = zeros(numInstances,numMethods);
totalTime = zeros(numInstances,numMethods);
status = strings(numInstances,numMethods);

for i = 1:numInstances
    instanceNames(i) = string(BIPtable{i}.pinfo.name);
    numVariables(i) = BIPtable{i}.pinfo.n;
    originalOrder(i) = numVariables(i) + 1;
    for j = 1:numMethods
        result = BIPtable{i}.FRA{j};
        reducedOrder(i,j) = result.r;
        orderRatio(i,j) = result.r/originalOrder(i);
        auxiliaryLPTime(i,j) = result.fra_time;
        basisTime(i,j) = result.frv_time;
        totalTime(i,j) = result.total_time;
        status(i,j) = string(result.status);
    end
end

reduced = reducedOrder < originalOrder;
successful = status == "OPTIMAL" | status == "SUBOPTIMAL";
timeLimit = status == "TIME_LIMIT";
otherFailure = ~successful & ~timeLimit;

reducedCount = sum(reduced,1)';
reducedPercent = 100*reducedCount/numInstances;
meanOrderRatioAll = mean(orderRatio,1)';
medianOrderRatioReduced = NaN(numMethods,1);
for j = 1:numMethods
    if any(reduced(:,j))
        medianOrderRatioReduced(j) = median(orderRatio(reduced(:,j),j));
    end
end
timeLimitCount = sum(timeLimit,1)';
otherFailureCount = sum(otherFailure,1)';

reductionSummary = table(methodNames,reducedCount,reducedPercent, ...
    meanOrderRatioAll,medianOrderRatioReduced,timeLimitCount, ...
    otherFailureCount,'VariableNames', ...
    {'Method','ReducedInstances','ReducedPercent','MeanOrderRatioAll', ...
    'MedianOrderRatioReduced','TimeLimits','OtherFailures'});

binLower = [0;2000;4000;6000;8000];
binUpper = [2000;4000;6000;8000;10000];
binLabels = ["n <= 2000"; "2000 < n <= 4000"; ...
    "4000 < n <= 6000"; "6000 < n <= 8000"; ...
    "8000 < n <= 10000"];
latexBinLabels = {'$n\leq 2000$', '$2000<n\leq 4000$', ...
    '$4000<n\leq 6000$', '$6000<n\leq 8000$', ...
    '$8000<n\leq 10000$'};
numBins = numel(binLower);
binCount = zeros(numBins,1);
medianTime = NaN(numBins,numMethods);
meanTime = NaN(numBins,numMethods);

for b = 1:numBins
    if b == 1
        inBin = numVariables <= binUpper(b);
    else
        inBin = numVariables > binLower(b) & numVariables <= binUpper(b);
    end
    binCount(b) = sum(inBin);
    if any(inBin)
        medianTime(b,:) = median(totalTime(inBin,:),1);
        meanTime(b,:) = mean(totalTime(inBin,:),1);
    end
end

timingSummary = table(binLabels,binCount, ...
    medianTime(:,1),meanTime(:,1),medianTime(:,2),meanTime(:,2), ...
    medianTime(:,3),meanTime(:,3),'VariableNames', ...
    {'VariableRange','Instances','AffineFRMedian','AffineFRMean', ...
    'PartialFRDMedian','PartialFRDMean','PartialFRDDMedian', ...
    'PartialFRDDMean'});

numRows = numInstances*numMethods;
Instance = strings(numRows,1);
NumVariables = zeros(numRows,1);
OriginalOrder = zeros(numRows,1);
Method = strings(numRows,1);
Status = strings(numRows,1);
ReducedOrder = zeros(numRows,1);
OrderRatio = zeros(numRows,1);
Reduced = false(numRows,1);
AuxiliaryLPTime = zeros(numRows,1);
BasisTime = zeros(numRows,1);
TotalTime = zeros(numRows,1);

row = 0;
for i = 1:numInstances
    for j = 1:numMethods
        row = row + 1;
        Instance(row) = instanceNames(i);
        NumVariables(row) = numVariables(i);
        OriginalOrder(row) = originalOrder(i);
        Method(row) = methodNames(j);
        Status(row) = status(i,j);
        ReducedOrder(row) = reducedOrder(i,j);
        OrderRatio(row) = orderRatio(i,j);
        Reduced(row) = reduced(i,j);
        AuxiliaryLPTime(row) = auxiliaryLPTime(i,j);
        BasisTime(row) = basisTime(i,j);
        TotalTime(row) = totalTime(i,j);
    end
end

instanceResults = table(Instance,NumVariables,OriginalOrder,Method,Status, ...
    ReducedOrder,OrderRatio,Reduced,AuxiliaryLPTime,BasisTime,TotalTime);

prefix = "miplib_" + string(timestamp);
files.reductionCSV = fullfile(resultsDir,prefix + "_reduction_summary.csv");
files.timingCSV = fullfile(resultsDir,prefix + "_timing_summary.csv");
files.instanceCSV = fullfile(resultsDir,prefix + "_instance_results.csv");
files.summaryMAT = fullfile(resultsDir,prefix + "_summary.mat");
files.latexFile = fullfile(resultsDir,prefix + "_tables.tex");

writetable(reductionSummary,files.reductionCSV);
writetable(timingSummary,files.timingCSV);
writetable(instanceResults,files.instanceCSV);
save(files.summaryMAT,'reductionSummary','timingSummary', ...
    'instanceResults','-v7.3');

write_latex_tables(files.latexFile,reductionSummary,timingSummary, ...
    latexMethodNames,latexBinLabels);

end

function write_latex_tables(filename,reductionSummary,timingSummary, ...
    latexMethodNames,latexBinLabels)

fid = fopen(filename,'w');
if fid < 0
    error('Unable to open LaTeX output file: %s',filename);
end
fileCleanup = onCleanup(@() fclose(fid));
rowEnd = '\\';

fprintf(fid,'%s\n','% Generated by generate_miplib_tables.m.');
fprintf(fid,'%s\n','% Each timing entry is median (mean) total preprocessing time in seconds.');
fprintf(fid,'%s\n','');

fprintf(fid,'%s\n','\begin{table}[htbp]');
fprintf(fid,'%s\n','\centering');
fprintf(fid,'%s\n','\small');
fprintf(fid,'%s\n','\begin{tabular}{@{}lrrrrr@{}}');
fprintf(fid,'%s\n','\toprule');
fprintf(fid,'%s & %s & %s & %s & %s & %s %s\n', ...
    'Method','\shortstack{Reduced\\instances}', ...
    '\shortstack{Mean ratio\\(all)}', ...
    '\shortstack{Median ratio\\(reduced)}', ...
    '\shortstack{Time\\limits}', ...
    '\shortstack{Other\\failures}',rowEnd);
fprintf(fid,'%s\n','\midrule');
for j = 1:height(reductionSummary)
    reducedCell = sprintf('%d (%.1f%%)', ...
        reductionSummary.ReducedInstances(j), ...
        reductionSummary.ReducedPercent(j));
    reducedCell = strrep(reducedCell,'%','\%');
    fprintf(fid,'%s & %s & %s & %s & %d & %d %s\n', ...
        latexMethodNames{j},reducedCell, ...
        format_ratio(reductionSummary.MeanOrderRatioAll(j)), ...
        format_ratio(reductionSummary.MedianOrderRatioReduced(j)), ...
        reductionSummary.TimeLimits(j), ...
        reductionSummary.OtherFailures(j),rowEnd);
end
fprintf(fid,'%s\n','\bottomrule');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s%d%s\n', ...
    '\caption{Matrix-order reductions across the ', ...
    sum(timingSummary.Instances), ...
    [' MIPLIB mixed-binary instances. Ratios are reduced order ', ...
    'divided by original order.}']);
fprintf(fid,'%s\n','\label{froverview}');
fprintf(fid,'%s\n','\end{table}');
fprintf(fid,'%s\n','');

fprintf(fid,'%s\n','\begin{table}[htbp]');
fprintf(fid,'%s\n','\centering');
fprintf(fid,'%s\n','\small');
fprintf(fid,'%s\n','\begin{tabular}{lrrrr}');
fprintf(fid,'%s\n','\toprule');
fprintf(fid,'%s & %s & %s & %s & %s %s\n', ...
    'Number of variables','Instances','Affine FR', ...
    'Partial FR ($\mathcal{D}$)','Partial FR ($\mathcal{DD}$)',rowEnd);
fprintf(fid,'%s\n','\midrule');
for b = 1:height(timingSummary)
    fprintf(fid,'%s & %d & %s & %s & %s %s\n', ...
        latexBinLabels{b},timingSummary.Instances(b), ...
        format_time(timingSummary.AffineFRMedian(b), ...
            timingSummary.AffineFRMean(b)), ...
        format_time(timingSummary.PartialFRDMedian(b), ...
            timingSummary.PartialFRDMean(b)), ...
        format_time(timingSummary.PartialFRDDMedian(b), ...
            timingSummary.PartialFRDDMean(b)),rowEnd);
end
fprintf(fid,'%s\n','\bottomrule');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n',['\caption{Preprocessing time in seconds by problem size. ', ...
    'Each entry is median (mean) total preprocessing time.}']);
fprintf(fid,'%s\n','\label{frtime}');
fprintf(fid,'%s\n','\end{table}');

end

function value = format_ratio(x)
if isnan(x)
    value = '--';
else
    value = sprintf('%.4f',x);
end
end

function value = format_time(medianValue,meanValue)
if isnan(medianValue) || isnan(meanValue)
    value = '--';
else
    value = sprintf('%.2f (%.2f)',medianValue,meanValue);
end
end
