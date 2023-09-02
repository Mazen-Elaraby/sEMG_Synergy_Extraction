function [rec_data] = sEMG_reconstruction(synergy_mat,time_act)

[NoC, ~, ~] = size(synergy_mat);
[~ , num_samples, L] = size(time_act);

rec_data = zeros(num_samples, NoC, L, 'single');

for i=1:L
    rec_data(:,:,i) = (synergy_mat(:,:,i) * time_act(:,:,i))';
end

end