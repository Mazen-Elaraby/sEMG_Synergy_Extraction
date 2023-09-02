function [MVC_normalized] = sEMG_preprocessing(raw_singal)

Fs = 2000; %Sampling Freq 2 KHz

[num_samples, NoC] = size(raw_singal); %No. of Channels corr. to No. of muscles

%Band-pass filtering

Fn = Fs/2; %Nyquist Freq
LP = 25; %low cut-off
HP = 350; %high cut-off
[b, a] = butter(4,[LP,HP]/Fn,'bandpass');
emg_filt = zeros(num_samples, NoC,'single'); %band-pass filtered data
for i=1:NoC
 emg_filt(:,i) = filtfilt(b,a,raw_singal(:,i));
end

%Full-wave Rectification

emg_rect = abs(emg_filt);

%RMS Enevelope Generation

emg_envelope = zeros(num_samples,NoC,'single');
window = 50;
for i=1:NoC
 emg_envelope(:,i) = sqrt(movmean(emg_rect(:,i).^2,window));
end

%MVC Normalization

MVCs = max(abs(raw_singal));
MVC_normalized = zeros(num_samples, NoC,'single');
for i=1:NoC
 MVC_normalized(:,i) = emg_envelope(:,i) ./ MVCs(1,i);
end

end