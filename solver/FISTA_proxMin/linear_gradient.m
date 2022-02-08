function [f, g] = linear_gradient(x,A,Ah,b)
u = A(x); clear x
r = u-b; clear u
g = Ah(r);
f = norm(r,'fro');