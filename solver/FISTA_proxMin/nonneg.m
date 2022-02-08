function [denoised, norm_out] = nonneg(x)
% denoised = max(x,0);
norm_out = 0;

cc = (size(x,2)/4+1):(3*size(x,2)/4);
rc = (size(x,1)/4+1):(3*size(x,1)/4);
x_temp = x(rc,cc);
denoised = zeros('like',x);
denoised(rc,cc) = x_temp;
denoised = max(x,0);


