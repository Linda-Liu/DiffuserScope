function Atb = A_adj_lensless_3d(h,x,gputrue)
% Atb = A_adj_lensless_3d(h,x,crop,pad,gputrue)
% inputs: see A
%
if gputrue
    Atb = gpuArray(zeros(size(h)));
else
    Atb = zeros(size(h));
end
Bp = fft2(x);
%parfor m = 1:size(h,3)
for m = 1:size(h,3)
    H = conj(fft2(h(:,:,m)));
    Atb(:,:,m) = ifftshift(real(ifft2(H.*Bp)));
end