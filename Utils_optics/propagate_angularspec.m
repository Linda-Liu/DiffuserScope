function [uout] = propagate_angularspec(uin,L,lambda,z)
% PROPAGATE_ANGULARSPEC   
% Propagate a wavefield using angular spectrum(Rayleigh-Sommerfeld diffraction)
%
% Parameters:
%   u1 - source plane field
%   L - source and observation plane side length
%   lambda - wavelength
%   z - propagation distance(meters)
%   u2 - observation plane field

  [M,N]=size(uin);
  dx=L/M;
  k=2*pi/lambda;
  
  U1 = fft2(uin);
  clear uin
  
  % Compute the frequency sampling rate and frequency bin size.
  fs = 1/dx;
  df = fs/M;

  
  if mod(M,2) == 0
    fx = [ [0:(M/2-1)]*df fs/2 [-(M/2-1):-1]*df ];  % Frequencies for FFT with even number of samples
  else
    fx = [ (0:(M/2-0.5))*df (-(M/2-0.5):-1)*df ];  % Frequencies for FFT with odd number of samples
  end
  
  [FX,FY]=meshgrid(fx,fx);
  H=exp(1j*k*z*sqrt(1-lambda^2*(FX.^2+FY.^2)));  %transfer function
  if strcmp(class(U1), 'gpuArray')
      H = gpuArray(H);
  end
  
  U2 = H.*U1;
  clear U1 H
  uout = ifft2(U2);            %inverse fft, center obs field
  clear U2
%  imagesc(abs(uout));axis image;title(z);
end
