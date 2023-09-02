function [Perf_Anal] = sEMG_perf_anal(org_data,rec_data)

%Root-mean-square error (RMSE)

RMSE_rec = norm(org_data - rec_data,'fro')/sqrt(size(org_data,1)*size(org_data,2));

%Variance Accounted For

avg_vaf = mean(vaf(org_data,rec_data));

%Subspace

theta = rad2deg(subspace(org_data,rec_data));

%Stationary Hjorth Parameters-Time-Domain based Power Spectrum Descriptors

[org_Activity, org_Mobility, org_Complexity] =hjorth(org_data);
[rec_Activity, rec_Mobility, rec_Complexity] =hjorth(rec_data);

RMSE_Activity = norm(org_Activity - rec_Activity,'fro')/sqrt(length(org_Activity));
RMSE_Mobility = norm(org_Mobility - rec_Mobility,'fro')/sqrt(length(org_Mobility));
RMSE_Complexity = norm(org_Complexity - rec_Complexity,'fro')/sqrt(length(org_Complexity));

RMSE_Hjorth = [RMSE_Activity, RMSE_Mobility, RMSE_Complexity];

%KL Divergence

KL = KL_Divergence(org_data,rec_data);

%Performance Analysis Data:

Perf_Anal = [RMSE_rec, avg_vaf, theta, RMSE_Hjorth, KL];

end