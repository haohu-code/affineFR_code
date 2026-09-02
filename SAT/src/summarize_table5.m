function summary = summarize_table5(EXPoutputs,folderNames)
% SUMMARIZE_TABLE5  Aggregate rank-only affine-FR results for Table 5.

arguments
    EXPoutputs cell
    folderNames cell
end

tableOrder = table5_families();

if numel(EXPoutputs) ~= numel(folderNames)
    error('The result and dataset-name arrays must have the same length.');
end

if ~isequal(sort(folderNames(:)), sort(tableOrder))
    error(['Table 5 requires exactly these ten SAT families: ', ...
        strjoin(tableOrder', ', '), '.']);
end

numFamilies = numel(tableOrder);
numInstances = zeros(numFamilies,1);
numReduced = zeros(numFamilies,1);
averageOriginalOrder = zeros(numFamilies,1);
averageReducedOrder = zeros(numFamilies,1);
averageRatio = zeros(numFamilies,1);
averageTimeSeconds = zeros(numFamilies,1);

for i = 1:numFamilies
    sourceIndex = find(strcmp(folderNames,tableOrder{i}),1);
    familyResults = EXPoutputs{sourceIndex};
    familyResults = familyResults(~cellfun(@isempty,familyResults));

    if isempty(familyResults)
        error('No eligible instances were processed for family %s.', ...
            tableOrder{i});
    end

    originalOrder = cellfun(@(x) x.info.n + 1, familyResults);
    reducedOrder = cellfun(@(x) x.affineFRA.r, familyResults);
    preprocessingTime = cellfun(@(x) x.affineFRA.total_time, familyResults);

    numInstances(i) = numel(familyResults);
    numReduced(i) = sum(reducedOrder < originalOrder);
    averageOriginalOrder(i) = mean(originalOrder);
    averageReducedOrder(i) = mean(reducedOrder);
    averageRatio(i) = mean(reducedOrder ./ originalOrder);
    averageTimeSeconds(i) = mean(preprocessingTime);
end

summary = table(string(tableOrder), numInstances, numReduced, ...
    averageOriginalOrder, averageReducedOrder, averageRatio, ...
    averageTimeSeconds, 'VariableNames', {'Family', 'NumInstances', ...
    'NumReduced', 'AverageOriginalOrder', 'AverageReducedOrder', ...
    'AverageRatio', 'AverageTimeSeconds'});

end
