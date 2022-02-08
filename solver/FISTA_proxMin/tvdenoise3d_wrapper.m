function [denoised, norm_out] = tvdenoise3d_wrapper(x,tau,lambda2,niters,minval,maxval)
denoised = tvdenoise(x,2/tau,niters,lambda2);
denoised = min(max(denoised,minval),maxval);
norm_out = tau*TVnorm3d(denoised,lambda2);