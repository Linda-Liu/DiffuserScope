function [Xguess,loss] = deconvRLnew(forwardFUN, backwardFUN, Htf, maxIter, Xinit, b, options)
%  forwardFUN(Xguess) = H Xguess
%  backwardFUN(b)= Ht b
%  Xinit -  Initial guess
%  b      -  Measurement
%  Linda Liu 08/2019
if ~isfield(options,'xsize')
    options.xsize = size(Xinit);
end
if ~isfield(options,'save_progress')
    options.save_progress = 0;
end
if ~isfield(options,'color_map')
    options.color_map = parula;
end
if options.save_progress
    if ~isfield(options,'progress_file')
        options.progress_file = [options.savePath 'RL_progress.avi'];
    end
    if exist(options.progress_file,'file')
        overwrite_mov = input('video file exists. Overwrite? y to overwrite, n to abort.');
        if strcmpi(overwrite_mov,'n')
            new_mov_name = input('input new name (no extension): ');
            options.progress_file = [new_mov_name,'.avi'];
        end
    end
    options.vidObj = VideoWriter(options.progress_file);
    options.vidObj.Quality = 100;
    options.vidObj.FrameRate = 1;
    open(options.vidObj);
end

options.fighandle = figure(1);clf

onevec = ones(size(b));
normalization = backwardFUN(onevec); clear onevec
if min(normalization(:))<0
    normalization = normalization +0.1;
end
b_norm = norm(b,'fro');
nvoxels = numel(Xinit);

Xguess = Xinit; clear Xinit
for step_num=1:maxIter
    tic;
%     Xguess_prev = Xguess;
    
    %%%%%%%Original RL update%%%%%
    %{
    HXguess = forwardFUN(Xguess);
    HXguessBack = backwardFUN(HXguess);
    errorBack = Htf./HXguessBack;
    Xguess = Xguess.*errorBack; 
    %}
    
    % RL Update multiplication
    % https://github.com/scikit-image/scikit-image/blob/master/skimage/restoration/deconvolution.py
    %{
    clear Htf
    HXguess = forwardFUN(Xguess);
    relativs_blur= b./HXguess;
    residual_norm = norm(b-HXguess,'fro') / b_norm; clear HXguess
%     relativs_blur(HXguess<prctile(HXguess,1,'all')) = 0;
    Xguess = Xguess.*backwardFUN(relativs_blur)./normalization;
    Xguess(find(isnan(Xguess))) = 0;
    clear relativs_blur
    %}
    
    % RL Update additive
    %
    HXguess = forwardFUN(Xguess);
    relativs_blur= b./HXguess;
    Xguess = Xguess+ 1e-4*(backwardFUN(relativs_blur)-normalization);
    %}
    
    % Compute residual
%     residual_norm = norm(b-HXguess,'fro') / b_norm;
%     X_diff = reshape(Xguess-Xguess_prev,[nvoxels,1]);
%     Update_norm = norm(X_diff) / nvoxels;
    
    ttime = toc;
%     disp(['  iter ' num2str(step_num) ' | ' num2str(maxIter) ', took ' num2str(ttime) ' secs, Residual Norm ' num2str(residual_norm) ', Update Norm ' num2str(Update_norm)]);
    disp(['  iter ' num2str(step_num) ' | ' num2str(maxIter) ', took ' num2str(ttime) ' secs, Residual Norm ' num2str(residual_norm)]);
    loss(step_num) = residual_norm;
    
    if ~mod(step_num,options.disp_fig_interval)
        if options.disp_figs
            draw_figures(Xguess,options);
            %suptitle(['iter ', num2str(step_num)]);
            figure(2)
            plot(loss(2:end))
            ylabel('loss');xlim([1 maxIter])
            drawnow
        end
        if options.save_progress
            frame = getframe(options.fighandle);
            writeVideo(options.vidObj,frame);
        end
    end
end
if options.save_progress
    close(options.vidObj);
end

if isfield(options,'psf_norm')
    Xguess = re_normalize(Xguess,options);
end
return

function xout = re_normalize(xk,options)
        xout = bsxfun(@times,xk,reshape(1./options.psf_norm,1,1,[]));
end

function draw_figures(xk,options)
    if isfield(options,'psf_norm')
        xk = re_normalize(xk,options);
    end
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
    im1 = squeeze(max(xk,[],3));
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
end