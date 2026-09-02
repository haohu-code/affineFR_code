function [weights,record] = load_table6_objective(instanceName)
%LOAD_TABLE6_OBJECTIVE Load the fixed linear objective used in Table 6.

name = regexprep(char(string(instanceName)),'\.cnf$','');

persistent table6Objectives
if isempty(table6Objectives)
    srcDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(srcDir);
    dataFile = fullfile(projectRoot,'data','table6_objectives.mat');
    savedData = load(dataFile,'table6Objectives');
    table6Objectives = savedData.table6Objectives;
end

index = find(strcmp({table6Objectives.instance},name),1);
if isempty(index)
    error('No fixed Table 6 objective is stored for instance %s.',name);
end

record = table6Objectives(index);
weights = record.weights(:);

end
