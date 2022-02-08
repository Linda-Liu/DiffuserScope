function [denoised, norm_out] = tvdenoise2d_wrapper(x,tau,niters,minval,maxval)
denoised = tvdenoise(x,2/tau,niters);

denoised = min(max(denoised,minval),maxval);
norm_out = tau*TVnorm(denoised);