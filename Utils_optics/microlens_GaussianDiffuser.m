function [phase] = microlens_GaussianDiffuser(M,L,lambda,D,ns,zf)
% MICROLENS_GAUSSIANDIFFUSER generate gaussian diffuser transfer function
%
% Parameters: 
%   M - input field size
%   L - side length of the diffuser(meters)
%   lambda - wavelength(meter)
%   D - lenslet pitch on average (meters)
%   ns - number of lenslets in each linear dimension
%   zf - focal length on average

%   degree - diffusing angle(diffuser deflection angle)
%  "focal length"=D^2*delta_n/(degree*pi/180) (meters)

% Save booleans
save_m = false;

delta_n = .5; % refractive index difference
dx=L/M; %simulation pixel size

%% generate real-space Gaussuian filter
sigma = D/6; % standard deviation of Gaussian(meters)
A = sigma^2/delta_n/zf*2.4;
fprintf('A is %0.4e meter\n', A)

hsize = round(min(sigma/dx*12,M)); % gaussian filter size to be 5 sigma
K = fspecial('gaussian',hsize,sigma/dx); %fspecial('gaussian',HSIZE,SIGMA)

%% apply gaussian filter to a random surface
seeds = rand(M);  %uniformly distrubuted in the interval (0,1)
surface_height = ifft2(fft2(K,M,M).*(fft2(seeds))); %unit meter
surface_height = surface_height-min(surface_height(:));
%adjust the height
surface_height = surface_height/prctile(surface_height(:),99.99) * A;
% figure;
imagesc(surface_height*delta_n*2*pi/lambda);axis image;

phase = exp(1j*surface_height*delta_n*2*pi/lambda);
% figure;imagesc(angle(phase));axis image;

if save_m
    fname = '../diffuser_phase';
    save(strcat(fname, '.m'),'phase');
end

end