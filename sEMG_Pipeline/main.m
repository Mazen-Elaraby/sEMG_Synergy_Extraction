%% Collecting Data
data_paths = dir('**/*E1*.mat');

num_subjs = length(data_paths);

emg_data = cell([1,num_subjs]);
for i=1:num_subjs
    emg_data{i} = load(fullfile(data_paths(i).folder, data_paths(i).name)).emg;
end
%% sEMG Data Preprocessing

preprocessed_data = cell([1,num_subjs]);
for i=1:num_subjs
    preprocessed_data{i} = sEMG_preprocessing(emg_data{i});
end
%% Synergy Extraction

addpath("FastICA_25","Helper_Functions");
no_syn = 4;

synergy_mat = cell([1,num_subjs]);
time_act = cell([1,num_subjs]);

for i=1:num_subjs
   [synergy_mat{i},time_act{i}] = sEMG_synergy_extraction(preprocessed_data{i},no_syn);
end
%% sEMG Data Reconstruction 

rec_data = cell([1,num_subjs]);

for i=1:num_subjs
    rec_data{i} = sEMG_reconstruction(synergy_mat{i}, time_act{i});
end

%% Performance Analysis

%Performace Analysis is based on 7 metric and tested on 4 synergy
%extraction methods
num_metrics = 7;
num_methods = 4;

perf_anal = zeros(num_methods, num_metrics, num_subjs,'single');

for i=1:num_subjs
    for j=1:num_methods
        perf_anal(j,:,i) = sEMG_perf_anal(preprocessed_data{i}, rec_data{i}(:,:,j));
    end
end

avg_perf_anal = mean(perf_anal,3);
%% Save to mat file

filename = "sEMG_Processing_Pipeline_Results";
save(filename);
save(filename,"synergy_mat","time_act","avg_perf_anal");


