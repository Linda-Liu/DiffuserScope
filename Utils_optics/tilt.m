function[uout]=tilt(uin,L,lambda,alpha,theta, usegpu)
% tilt phasefront
% uin - input field
% L - side length (meter)
% lambda - wavelength (meter)
% alpha - tilt angle (radian)
% theta - rotation angle (x axis 0)
% Computational Fourier Optics- David Voelz

[M,N]=size(uin); %get input field array size
dx=L/M; %sample interval
k=2*pi/lambda; %wavenumber

if usegpu
    x=gpuArray(-L/2+dx/2:dx:L/2-dx/2); %coords
else
    x=(-L/2+dx/2:dx:L/2-dx/2); %coords
end
[X,Y]=meshgrid(x,x);

uout=uin.*exp(1i*k*(X*cos(theta)+Y*sin(theta))*tan(alpha)); %apply tilt
clear x X Y uin
% imagesc(angle(uout));axis image
end
