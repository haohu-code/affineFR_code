function clean_data()
% select only mixed binary optimization problems

% get the paths of the MIPLIB directory and its data
scriptDir = fileparts(mfilename('fullpath'));
miplibDir = fileparts(scriptDir);
dataDir = fullfile(miplibDir,'data');

% load the data
files = dir(dataDir);
prob_list = cell(1000,1);
k = 0;
for j = 1:length(files)
    if (files(j).isdir == true)
        continue
    end
    % get the extension
    dot_idx = find(files(j).name == '.');

    s = files(j).name(dot_idx+1:end);
    if strcmp(s,'mps')
        k = k + 1;
        prob_list{k} = files(j);
    end
end
prob_list = prob_list(1:k);

% sort the instances
psize = zeros(k,1);
for i = 1:k
    %
    p_path = [prob_list{i}.folder '/' prob_list{i}.name];
    prob = gurobi_read(p_path);
    psize(i) = length(prob.obj);
    prob_list{i}.n = psize(i);
%     keyboard
end
[~,idx] = sort(psize);
prob_list = prob_list(idx);

for i = 1:k
    prob_list{i} = rmfield(prob_list{i},'folder');
    prob_list{i} = rmfield(prob_list{i},'date');
    prob_list{i} = rmfield(prob_list{i},'bytes');
    prob_list{i} = rmfield(prob_list{i},'isdir');
    prob_list{i} = rmfield(prob_list{i},'datenum');
end
save(fullfile(miplibDir,'prob_list.mat'),'prob_list')


end
