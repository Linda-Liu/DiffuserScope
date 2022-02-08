% -------------------------------------------------------------------
% Reconstruc a video from Fourier Light Field microscope
%
% -------------------------------------------------------------------
addpath('/solver/')
addpath('/ForwardModel/')

%% propare PSF stack

indir = 'D:\Linda\200909_1um_bead_hydra\1um_calib_bead\psf_1\psf-bg';
psf_step = 4;
psf_start = 1;
psf_end = 400;

if exist(indir, 'dir')
    infiles_struct = dir(fullfile(indir, '/*.tif*'));
    [~, order] = sort({infiles_struct(:).name});
    infiles_struct = infiles_struct(order);
else
    disp('indir does not exist');
    return
end

infiles_struct = infiles_struct(psf_start:psf_step:psf_end);
h = zeros(2048,2048,length(infiles_struct));
for i=1:length(infiles_struct)
    hh = double(imread(fullfile(indir, infiles_struct(i).name), 'tiff'));
    h(:,:,i)=hh;
end

%crop center
center = [381,1800,381,1800];
h = h(center(1):center(2),center(3):center(4),:);

% downsample
lateral_downsample = 4;
for n = 1:log2(lateral_downsample)
    h = 1/4*(h(1:2:end,1:2:end,:)+h(1:2:end,2:2:end,:) + ...
        h(2:2:end,1:2:end,:) + h(2:2:end,2:2:end,:));
end

% normalize
h_n = zeros(1,size(h,3)); %record norm of each slice, multiply each slice in the end

for n = 1:size(h,3)
    h_n(n) = norm(h(:,:,n),'fro');
%     h(:,:,n) = h(:,:,n)/h_n(n);
end

h = h/max(h_n);

%% prepare recon frames
indir = 'D:\Linda\200909_1um_bead_hydra\hydra\z-100_10\frames';
img_step = 1;
img_start = 1;
img_end = 200;

if exist(indir, 'dir')
    infiles_struct = dir(fullfile(indir, '/*.tif*'));
    [~, order] = sort({infiles_struct(:).name});
    infiles_struct = infiles_struct(order);
else
    disp('indir does not exist');
    return
end

infiles_struct = infiles_struct(img_start:img_step:img_end);

%% for loop: recon every frame
for i=1:length(infiles_struct)
% for i=1:size(L,1)
% fprintf('Reconstruct frame %d \n',i+img_start-1)
frame = i;
fprintf('Reconstruct frame %d \n',frame)
b = double(imread(fullfile(indir, infiles_struct(i).name), 'tiff'));
% b = squeeze(L(i,:,:)); %recon sparse/Low-rank component

% crop center
b = b(center(1):center(2),center(3):center(4));

% downsample
for n = 1:log2(lateral_downsample)
    b = 1/4*(b(1:2:end,1:2:end)+b(1:2:end,2:2:end) + ...
        b(2:2:end,1:2:end) + b(2:2:end,2:2:end));
end

%normalize
norm_b = max(b(:));
b = b/norm_b;

usegpu = 1;
f=figure(1);clf
M = 710;
dx = 1*lateral_downsample; %in um
mag = 6.5;
dz = 1 * psf_step; %in um

% Deconvolution using PSFs 
options.disp_crop = @(x)gather(x(floor(size(x,1)*0.3):floor(size(x,1)*0.8),...
    floor(size(x,2)*0.25):floor(size(x,2)*0.75),5:end));    %Function to be applied to data before displaying    
%options.disp_crop = @(x)x;
options.disp_fig_interval = 1;%options.maxIter;
options.fighandle = f;
options.disp_figs = 0;
options.disp_gamma = 1;
options.color_map = 'parula';
options.save_progress = 0; %record iterations in video
options.save_every = 0; % Save image(.fig) and stack(.mat) every N iterations
% options.savePath = 'D:\Linda\210701_video_NMF_SLdecompose\210701_video_NMF_SLdecompose\SLdecompose\200909_hydra_z-100_10\L_lambda0p1_recon\'; %path to save .fig and .avi
options.savePath = 'D:\Linda\200909_1um_bead_hydra\Reconstruction\z-100_10\';
options.xsize = size(h);
options.known_input = 0;   % If input is known, compute PSNR
%options.psf_norm = h_n;
options.tv_tau = 1e-3;
options.lambda2 = dx/dz;

A3d = @(x)A_lensless_3d(h, x, usegpu);
Aadj_3d = @(x)A_adj_lensless_3d(h, x, usegpu);

%%
%------------------- RL Reconstruction-------------------------------------
%
options.maxIter = 200;  %Number of iterations 125 for ori, 80 for S component
if usegpu    
    Htf = Aadj_3d(gpuArray(b));
else
    Htf = Aadj_3d(b);
end
Htf(Htf<0) = 0;
% Xguess = Htf
Xhat = deconvRLTVnew(A3d, Aadj_3d, Htf, options.maxIter, Htf, b, options);
Xhat = Xhat*norm_b; %restore normalization
Xhat = gather(Xhat);

% save 
input_name = ['frame', num2str(frame)];
output_name = [options.savePath,'xhat/',input_name,'_Xhat','.mat'];
save(output_name, 'Xhat','-v7.3');

output_img = [options.savePath,'tiff/',input_name,'_Xhat','.tif'];
saveas(options.fighandle, output_img);
end
%}

%%
%------------------- FISTA Reconstruction----------------------------------
%
% deconvolution options
options.maxIter = 250; 
options.stepsize = 5e-6;
options.convTol = 1e-30;   %stop if norm(f(xk+1)-f(xk)) is small
options.residTol = 1e-10;   %Stopping tolerance norm(Ax-b)
options.momentum = 'nesterov';  %linear or nesterov. There's pretty much no reason not to use nesterov
options.tv_tau = 1e-11;
% used for anisotropic TV, lambda2=pixel_size_xy/pixel_size_z  
tv_lambda2 = dx/dz; % set to 1 if isotropic
tau = 1e-13; %for natural sparse
% Define gradient handle
GradErrHandle = @(x) linear_gradient(x,A3d,Aadj_3d,b);

% Prox handle: non-negativity, crop the center microlens for mla recon
% prox_handle = @(x)nonneg(x);
% prox_handle = @(x)soft3d(x,tau);
prox_handle = @(x)tvdenoise3d_wrapper(x,options.tv_tau,tv_lambda2,15,0,inf);

Xguess = Aadj_3d(b);
[Xhat, funvals] = proxMin(GradErrHandle,prox_handle,Xguess,b,options);
clear Xguess

Xhat = Xhat*norm_b; %restore normalization
Xhat = gather(Xhat);

draw_figures(Xhat_out,options)

%%
input_name = ['frame', num2str(frame)];
output_name = [options.savePath,'xhat/',input_name,'_Xhat','.mat'];
save(output_name, 'Xhat_out','-v7.3');

output_img = [options.savePath,'tiff/',input_name,'_Xhat','.tif'];
saveas(options.fighandle, output_img);

end
%%
function draw_figures(xk,options)
xk = options.disp_crop(xk);
set(0,'CurrentFigure',options.fighandle)
xk = max(xk,0);  % add non-negativity
xk = xk.^options.disp_gamma;

if numel(options.xsize)==2
    imagesc(options.disp_crop(xk))
    axis image
    colorbar
    colormap(options.color_map);
    %caxis(gather([prctile(xk(:),.1) prctile(xk(:),90)]))
elseif numel(options.xsize)==3
    xk = gather(xk);
    set(0,'CurrentFigure',options.fighandle)
    subplot(1,3,1)
    im1 = imrotate(squeeze(max(xk,[],3)),180);
    imagesc(im1);
    hold on
    axis image
    colormap(options.color_map);
    %colorbar
    caxis([0 prctile(im1(:),99.9)])
    title('xy');
    set(gca,'fontSize',8)
    axis off
    hold off
    set(0,'CurrentFigure',options.fighandle)
    
    subplot(1,3,2)
    im2 = squeeze(max(xk,[],1));
    imagesc(im2);
    hold on    
    %axis image
    colormap(options.color_map);
    %colorbar
    title('yz');
    set(gca,'fontSize',8)
    caxis([0 prctile(im2(:),99.9)])
    axis off
    hold off
    drawnow
    set(0,'CurrentFigure',options.fighandle)
    
    subplot(1,3,3)
    im3 = squeeze(max(xk,[],2));
    imagesc(im3);
    hold on
    %axis image
    colormap(options.color_map);
    %colorbar   
    title('xz');
    set(gca,'fontSize',8)
    caxis([0 prctile(im3(:),99.9)]);
    axis off
    hold off
    

elseif numel(options.xsize) == 4
    xkr = reshape(xk,options.xsize);
    subplot(2,2,1)
    imagesc(transpose(squeeze(xkr(end,ceil(options.xsize(2)/2),:,:))))
    hold on
    axis image
    colorbar
    colormap gray
    caxis([0 prctile(xkr(:),99)]);
    hold off
    
    subplot(2,2,2)
    imagesc(transpose(squeeze(xkr(1,ceil(options.xsize(2)/2),:,:))))
    hold on
    axis image
    colorbar
    colormap gray
    caxis([0 prctile(xkr(:),99)]);
    hold off
    
    subplot(2,2,3)
    imagesc(transpose(squeeze(xkr(ceil(options.xsize(2)/2),1,:,:))))
    hold on
    axis image
    colorbar
    colormap gray
    caxis([0 prctile(xkr(:),99)]);
    hold off
    
    subplot(2,2,4)
    imagesc(transpose(squeeze(xkr(ceil(options.xsize(2)/2),end,:,:))))
    hold on
    axis image
    colorbar
    colormap gray
    caxis([0 prctile(xkr(:),99)]);
    hold off
    
    
elseif numel(options.xsize)==1
    
    plot(xk)
end
drawnow
end
