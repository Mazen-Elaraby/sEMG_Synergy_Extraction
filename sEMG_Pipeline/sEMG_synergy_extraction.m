function [W,H] = sEMG_synergy_extraction(emg_data,no_syn)

%Principal Component Analysis

[coeff,score]=pca(emg_data);
W_PCA = coeff(:,1:no_syn);
H_PCA = score(:,1:no_syn)';

%Non-negative Matrix factorization

[W_NMF,H_NMF] = nnmf(emg_data',no_syn);

%Independent Component Analysis

[H_ICA,W_ICA, ~] = fastica(emg_data','lastEig',no_syn,'verbose','off');

%Autoencoder
%Training the neural network
autoenc = trainAutoencoder(emg_data',no_syn);
%Extracting Synergies & Time Activations
W_AE = autoenc.DecoderWeights;
H_AE = encode(autoenc, emg_data');

%Stacking results
NoC = size(W_AE,1);
num_samples = size(H_AE,2);

W = zeros(NoC, no_syn, 4,'single');

W(:,:,1) = W_PCA;
W(:,:,2) = W_NMF;
W(:,:,3) = W_ICA;
W(:,:,4) = W_AE;

H = zeros(no_syn, num_samples, 4, 'single');

H(:,:,1) = H_PCA;
H(:,:,2) = H_NMF;
H(:,:,3) = H_ICA;
H(:,:,4) = H_AE;

end