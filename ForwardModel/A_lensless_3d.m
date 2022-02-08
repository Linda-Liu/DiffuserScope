function b = A_lensless_3d(h,x,gputrue)
% b = A_lensless_3d(h,x,pad,crop,gputrue)
% Computes convolutional forward operator for 3D fluorescence imaging
% with diffuser mask.
% Inputs:
% h: impulse response stack
% x: variable to convolve with h
% pad: funcion handle to do zero padding for handling circular conv
% crop: function handle to crop final result back to sensor size
%
% Output: estimate of sensor data
%Initialize empty array in frequency domain
if gputrue
    B = gpuArray(zeros(size(h,1),size(h,2),size(h,3)));
else
    B = zeros(size(h,1),size(h,2),size(h,3));
end

%parfor m = 1:size(h,3)
for m = 1:size(h,3)
    B(:,:,m) = fft2(x(:,:,m)).*fft2(h(:,:,m));
end
B = sum(B,3);

b = ifftshift(real(ifft2(B)));