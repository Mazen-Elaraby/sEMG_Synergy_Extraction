%% Choosing Hardware Acceleration Solution

%Not enough host-side memory is available to allocate the workspace for
%MLDIVIDE.?????

%Conclusions:
%General: single percision arrays, in-place processing, interpolate or
%not?(trade off between speed and interpolation problems?)
%Parfor: acessing cell or tensor pages 
%        built-in support or multi-process? (if available)
%GPU:   cell arrays cannot be stored on gpu
%           dealing w/ cellarray complexities or direct tensor (also access time)
%           pagefun cannot be used with cell arrays, no cellfun
%       for built-in support, use it or pagefun?
%acceleration solution:

%Calculating speedup acheived by each method relative to case 0 (NMF)
no_syn = 6;
no_subjs = 10;
%no_subjs = length(emg_data);

%creating test data
%for CPU processing
%cell array
test_data_C = cell([1,no_subjs]);
for i=1:no_subjs
    test_data_C{i} = randn(2000000,12,'single');
end
%Tensor/multi-dimensional array
test_data_T = randn(2000000,12,no_subjs,'single');


%% Profiling...
%case 0: No acceleration

%case 0.1: No acceleration, acessing matrices from cell arrays
tic
W_0C = cell([1,no_subjs]);
H_0C = cell([1,no_subjs]);
for i=1:no_subjs
    [W_0C{i},H_0C{i}] = nnmf(test_data_C{i}',no_syn);
end
t_0C = toc;
disp(['The unaccelerated case (Cellarr):', num2str(t_0C),' seconds'])

%case 0.2: No acceleration, acessing pages of a multi-dimensional array
tic
W_0T = zeros(12,no_syn,no_subjs,'single');
H_0T = zeros(no_syn,2000000,no_subjs,'single');
for i=1:no_subjs
    [W_0T(:,:,i),H_0T(:,:,i)] = nnmf(test_data_T(:,:,i)',no_syn);
end
t_0T = toc;
disp(['The unaccelerated case (Tensor):', num2str(t_0T),' seconds'])

%% case 1: Parfor
parpool;
delete(gcp('nocreate'));
%case 1.1: Parfor, accessing matrices from cell arrays
tic
W_1C = cell([1,no_subjs]);
H_1C = cell([1,no_subjs]);
parfor i=1:no_subjs
    [W_1C{i},H_1C{i}] = nnmf(test_data_C{i}',no_syn);
end
t_1C = toc;
disp(['Using Parfor (Cellarr):', num2str(t_1C),' seconds'])

%case 1.2: Parfor, accessing pages of a multi-dimensional array
tic
W_1T = zeros(12,no_syn,no_subjs,'single');
H_1T = zeros(no_syn,2000000,no_subjs,'single');
for i=1:no_subjs
    [W_1T(:,:,i),H_1T(:,:,i)] = nnmf(test_data_T(:,:,i)',no_syn);
end
t_1T = toc;
disp(['Using Parfor (Tensor):', num2str(t_1T),' seconds'])


%% case 2: GPU Computing
no_syn = 6;
no_subjs = 10;
D = gpuDevice;
reset(D);

%creating test data for GPU processing
%cell array
G_test_data_C = cell([1,no_subjs]);
for i=1:no_subjs
    G_test_data_C{i} = randn(2000000,12,'single','gpuArray');
end

%case 2.1: gpuArray, accessing matrices from cell array, serial cellfun
wait(D)
tic
W_2C = cell([1,no_subjs]);
H_2C = cell([1,no_subjs]);
for i=1:no_subjs
    [W_2C{i},H_2C{i}] = nnmf(G_test_data_C{i}',no_syn);
end
wait(D)
t_2C = toc;
disp(['GPU Acceleration (Cellarr):', num2str(t_2C),' seconds'])
reset(D);

%% 
no_syn = 6;
no_subjs = 10;
D = gpuDevice;
reset(D);
gpuDevice().AvailableMemory
%creating test data for GPU processing
%Tensor/multi-dimensional array
G_test_data_T = randn(2000000,12,no_subjs,'single','gpuArray');
gpuDevice().AvailableMemory

%case 2.2: gpuArray, acessing pages of a multi-dimensional array, pagefun
wait(D)
tic
W_2T = zeros(12,no_syn,no_subjs,'single','gpuArray');
H_2T = zeros(no_syn,2000000,no_subjs,'single','gpuArray');
gpuDevice().AvailableMemory

for i=1:no_subjs
    [W_2T(:,:,i),H_2T(:,:,i)] = nnmf(G_test_data_T(:,:,i)',no_syn);
end
wait(D)
t_2T = toc;
disp(['GPU Acceleration (Tensor):', num2str(t_2T),' seconds'])
gpuDevice().AvailableMemory
reset(D);




%GPU consideration: for built-in support: pagefun, or regular?
%try and understand arrayfun and pagefun first

