function [y, norm_out]= soft3d(x,tau)

y = max(abs(x) - tau, 0);
y = y./(y+tau) .* x;

norm_out = 0;
for n=1:size(x,3)
    norm_out = norm_out+ tau*norm(x(:,:,n),1);
end

