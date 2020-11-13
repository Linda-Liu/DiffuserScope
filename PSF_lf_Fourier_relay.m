function Isensor = PSF_lf_Fourier_relay(x,y,z,recipe,display,phasemask)
% Isensor = PSF_lf_Fourier_relay(0,0,0,'optical_recipe_Fourier_Olympus',1,rmla);

% PSF_lf_Fourier_relay Generate PSF for a Fourier space light field microscope,
% by calculating the sperical wavefront at back pupil plane, demagnify by a relay system,
% and propagate to the sensor.
%
% x,y,z, - point source locations(meters)
% recipe - optical recipe file
 
fprintf('Calculating Fourier LFM PSF of x=%0.2e, y=%0.2e, z=%0.2e \n',x,y,z)
run(recipe)
usegpu = 0;
if ~exist('display','var')
      display = 1;
end
%% POINT SOURCE 
I = 1;  % point source intensity
% fprintf('Total intensity at the point source: %0.2e\n', I)
lambda = emission_wavelength;
zi = objective_focal_length.^2/(-z);

%% SIMULATION OPTIONS
downsample_at_sensor = 0;      % Choose whether to downsample the final intensity pattern to the sensor resolution

L = cam_size*demag;    % Simulation size in the back pupil space(meters). 

% To avoid aliasing issues, we run the simulation at or above the nyquist rate.
% The sim_sampling parameter is set automatically to ensure critical sampling.
nyquist_rate = lambda * ulens_fnumber/ 2 ; %diffraction limited resolution divided by 2
sim_sampling = ceil(cam_pixel / nyquist_rate); %every sim_sampling simulation pixels represent one camera pixel

% M- pixel numbers in simulation
M = sim_sampling * cam_nu ;    % Number of pixels in Simulation (sim_sampling ensures we are above the nyquist rate)
df = L/M;  % pupil plane pixel size (m)
fprintf('Simulation size is %d pixels. Sampling rate is %0.2e (m) in sensor space.\n', M, df/demag)

%% Generate wavefront at pupil for on-axis point

%wavefront at pupil
u1 = wavefront_spherical(sqrt(I), lambda, medium_index, L, M, 0, 0, zi, usegpu);

%% apodization function
%back aperture size in frequency space
back_aperture_diameter = 2 * objective_na / lambda; 
% convert pupil size(meters) to frequency space unit(1/m)
coeff = back_aperture_diameter/objective_pupil; 

if usegpu
    x_ticks=gpuArray((-M*df/2+df/2:df:M*df/2-df/2)*coeff);
else
    x_ticks=(-M*df/2+df/2:df:M*df/2-df/2)*coeff;
end
[X,Y]=meshgrid(x_ticks, x_ticks); clear x_ticks
R = sqrt(X.*X + Y.*Y); clear X Y
thetas = asin(complex(R*lambda/medium_index));
apodization = sqrt(cos(thetas)); clear thetas
apodization(R>(back_aperture_diameter/2))=0; clear R

u2 = u1.*apodization;
%clear thetas X Y R u1 apodization
%% normalize the intensity at pupil so that every z layer has the same total enerygy
I_u2 = abs(u2).^2;
I_frac =(1-objective_focal_length/sqrt(objective_focal_length^2+objective_pupil^2/4))/2; % energy fraction of the solid angle extended by the back pupil
u2_normalize = u2/sqrt(sum(I_u2(:)))*sqrt(I*I_frac); %clear u2 I_u2 


%% Apply tilt to off-axis point
alpha = atan(sqrt(x^2+y^2)/objective_focal_length); % tilt angle (radian)
theta = atan2(y,x);% rotation angle (x axis 0)
u2_tilt = tilt(u2_normalize,L,lambda,alpha,theta, usegpu);
%clear u2 u2_normalize

%% De-magnify the relayed pupil
L = L/demag;

%% Apply phase mask at pupil
u3 = u2_tilt.*phasemask;
clear u2_tilt 

%% Propagate to the sensor plane
u4 = propagate_angularspec(u3,L,lambda,cam_distance);
       
%%            
Isensor = abs(u4).^2; clear u4
% fprintf('Total intensity at the image plane: %0.2e\n', sum(Isensor(:)));

%% downsample
if (sim_sampling > 1 && downsample_at_sensor)
  Isensor = imresize(Isensor,1/sim_sampling,'box');
  sim_sampling = 1;
  M = sim_sampling * cam_nu;
end

Isensor = gather(Isensor);

%% Display
if display
    clf
    gamma_correction = 0.5;        % Gamma correction (applied before displaying the intensity images below)
    % figure;
    % colordef black;
    % set(gcf,'Color','black')
    lims = [-L/2+df/2,L/2-df/2]/1e-6; %meters
    imagesc(lims, lims, Isensor.^gamma_correction);
    axis image
    colormap colormap_fire
    title(['\color{white}2D Intensity pattern at the sensor for point source at (',...
        num2str(round(x/1e-6)),',',num2str(round(y/1e-6)),',',num2str(z/1e-6),') micron.'])
    xlabel('\color{white}x (microns)')
    ylabel('\color{white}y (microns)')

    drawnow
end

end