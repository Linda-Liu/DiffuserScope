function u_out = abbe_sine_apodization(u_in, L, x, y, back_aperture_diameter, lambda, medium_index)
% ABBE_SINE_APODIZATION    
%
% Apply the apodization function for the Abbe sine condition.
% Inside pupil, apply P(theta)=sqrt(theta); outside pupil, P(theta)=0
% Note: only use this at the Fourier plane!

% L - Simulation size 
% x,y - center of the circular mask

  [M,~]=size(u_in);

  dx = L/M;
  x_ticks=-L/2+dx/2:dx:L/2-dx/2; %freq coords
  [X,Y]=meshgrid(x_ticks - x, x_ticks - y);
  R = sqrt(X.*X + Y.*Y);
%  R(find(R>(back_aperture_diameter/2))) = 0;

  thetas = asin(R*lambda/medium_index);
  apodization = sqrt(cos(thetas));
  apodization(R>(back_aperture_diameter/2))=0;
  u_out = apodization .* u_in;
end