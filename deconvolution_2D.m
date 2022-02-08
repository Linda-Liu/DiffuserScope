% -------------------------------------------------------------------
% Reconstruc a single 2D image from DiffuserCam
%
% -------------------------------------------------------------------
addpath('../solver/')
addpath('../solver/FISTA_proxMin')
addpath('../../Linda_opticslib/')
psfpath = '/Users/Linda/Documents/Research/LFM/data/220201_pdms/PSFstack/25.tiff';
bpath = '/Users/Linda/Documents/Research/LFM/data/220201_pdms/patrick_psf25.tiff';

ds = 8;

%% Prepare for reconstruction
h = double(imread(psfpath));
psf_norm = norm(h,'fro');
h = h / psf_norm;
b = double(imread(bpath));
b = b / psf_norm;

%downsample
h = imresize(h,1/ds,'box'); %binning
b = imresize(b,1/ds,'box'); %binning

%define crop and pad operators to handle 2D fft convolution
pad = @(x)padarray(x,[size(h,1)/2,size(h,2)/2],0,'both'); %create a handle to an anonymous function
cc = (size(h,2)/2+1):(3*size(h,2)/2);
rc = (size(h,1)/2+1):(3*size(h,1)/2);
crop = @(x)x(rc,cc);


%%
%------------------- Wiener deconvolution-------------------------------------
Xguess = deconvwnr(pad(b),pad(h),0.001); %deconvwnr(I,psf,nsr)
% Xguess(Xguess<0) = 0;
imagesc(crop(Xguess));

%% RL
Xguess = deconvlucy(pad(b),pad(h),50); %deconvlucy(I,psf,iter,damper)
% Xguess(Xguess<0) = 0;
imagesc(crop(Xguess));

%% regularized filter
Xguess = deconvreg(pad(b),pad(h),0.05); %deconvreg(I,psf,np)
imagesc(crop(Xguess));

%% Prepare for iteration
H = fft2(pad(h));
H_conj = conj(H);
HH = H.*H_conj;

% Define function handle for forward A(x)
A2d = @(x)real(crop(ifftshift(ifft2(H.*fft2(x)))));
Aadj_2d = @(x)real(ifftshift(ifft2(H_conj.*fft2(pad(x)))));
% A2d = @(x)lensless2d_forward(H, x);
% Aadj_2d = @(x)lensless2d_backward(H_conj, x);

% deconvolution options
options.stepsize = real(1.8 / max(HH(:)));
options.convTol = 1e-30;   %stop if norm(f(xk+1)-f(xk)) is small
options.maxIter = 500;  %Number of iterations
options.residTol = 1e-10;   %Stopping tolerance norm(Ax-b)
options.momentum = 'nesterov';  %linear or nesterov. There's pretty much no reason not to use nesterov
options.disp_crop = @(x)x;    %Function to be applied to data before displaying    
options.disp_fig_interval = 1;
options.fighandle = figure(1);
options.disp_figs = 1;
options.disp_gamma = 1;
options.color_map = 'gray';%'colormap_fire';
options.save_progress = 0; %record iterations in video
options.save_every = 20; % Save image(.fig) and stack(.mat) every N iterations
options.savePath = 'progress/'; %path to save .fig and .avi
options.xsize = size(b);
options.known_input = 0;   % If input is known, compute PSNR
options.tv_tau = 5e-4;
%%
%
%------------------- RL Reconstruction-------------------------------------
Htf = Aadj_2d(b);
% Htf(Htf<0) = 0;
[Xguess,loss] = deconvRLnew(A2d, Aadj_2d, Htf, options.maxIter, Htf, b, options);

%%
%------------------- FISTA Reconstruction----------------------------------
prox_handle = @(x)nonneg(x);
GradErrHandle = @(x) linear_gradient(x,A2d,Aadj_2d,b);
%Xinit = Aadj_2d(b);
Xinit = zeros(size(h)*2);
Xguess = proxMin(GradErrHandle,prox_handle,Xinit,b,options);
%}