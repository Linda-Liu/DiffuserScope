function [y, norm_out]= soft(x,tau)

y = max(abs(x) - tau, 0);
y = y./(y+tau) .* x;
norm_out = tau*norm(x,1);

