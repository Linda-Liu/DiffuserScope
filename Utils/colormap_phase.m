% __BEGIN_LICENSE__
%
% LFSIM: Light Field Microscope Wave Optics Simulations
%
% Please do not redistribute without permission.  See the README for
% points of contact & for more information about this software.
%
% Copyright (C) 2010-2013 Stanford University.
% All rights reserved.
%
% __END_LICENSE__

function p = colorphase(m)
% COLORMAP_PHASE   Red-blue phase colormap
% 
% COLORPHASE(M) returns an M-by-3 matrix containing a "phase" colormap.
% COLORMAP_PHASE, by itself, is the same length as the current figure's
% colormap. If no figure exists, MATLAB creates one.
%
% To add this colormap as a default map, use 'addpath' with the 
% directory containing 'colormap_phase.m'.
%
% To reset the colormap of the current figure use 'colormap(colormap_phase)'.
%
% see also:  HSV, GRAY, HOT, COOL, BONE, COPPER, FLAG, PINK, COLORMAP,
% RGBPLOT.
%

if nargin < 1
    m = size(get(gcf,'colormap'),1); 
end

%You can replace this M x 3 matrix with any matrix whose values range
%between 0 and 1 to create a new colormap file.  Use copy / paste to create
%a matrix like the one below, you do not have to add these values
%manually.  To create a new colormap, change 'cmap_mat' to the desired
%matrix, rename the function *and* the m-file from 'fire' to your desired
%colormap name.
x=[0:0.01:1]';
hue = sign(2*x - 1)/3 + 1/3;
sat = ones(size(x)) * 0.5;
val = abs(2*x - 1);
cmap_mat=hsv2rgb([hue sat val ]);

%interpolate values
xin=linspace(0,1,m)';
xorg=linspace(0,1,size(cmap_mat,1));

p(:,1)=interp1(xorg,cmap_mat(:,1),xin,'linear');
p(:,2)=interp1(xorg,cmap_mat(:,2),xin,'linear');
p(:,3)=interp1(xorg,cmap_mat(:,3),xin,'linear');





