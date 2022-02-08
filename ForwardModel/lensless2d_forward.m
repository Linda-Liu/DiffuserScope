function b = lensless2d_forward(H, x)
% A2d = @(x)real(ifftshift(ifft2(H.*fft2(x))));
X = fft2(x); clear x
foo = H.*X; clear H X
foo = ifft2(foo);
b = real(ifftshift(foo));

end