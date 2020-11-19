emission_wavelength = 510e-9;   % Wavelength of light to use in the simulation (m)

%% objective info N20X-PFH - 20X Olympus XLUMPLFLN Objective, 1.00 NA, 2.0 mm WD 
objective_mag = 20;             % Objective magnification
objective_na = 1.0;             % NA 
medium_index = 1.33;               % Index of refraction
f_tubelens = 180e-3;            % Focal length of the tube lens (m) (Nikon - 200mm, Olympus - 180mm)
d_tubelens = 180e-3;            % Actual distance between the tube lens and the back aperture plane (m)
f_relay = 48e-3;
objective_focal_length = f_tubelens/objective_mag;
objective_pupil = 18e-3;
demag = f_tubelens/f_relay;

%% camera 
% cam_pixel = 6.5e-06; %camera pixel size(meter) real number is 6.5e-06
cam_size = objective_pupil/demag; %camera sensor size in one dimension(meter)
cam_nu = 4500; %camera pixel numbers in one dimension
cam_pixel = cam_size/cam_nu;
cam_distance = 15.60*1e-3; %propagation distance from the phasemask to the sensor

%% microlens info
%ulens_fnumber=20;               
ulens_ns = 5; % number of lenslets in each linear dimension
ulens_pitch = objective_pupil/demag/5;   % Microlens pitch (m)  NA=1/2F-number
%ulens_focal_length = ulens_pitch*ulens_fnumber;     % Microlens focal length (m)
ulens_focal_length = cam_distance;
ulens_fnumber = ulens_focal_length/ulens_pitch;
ulens_profile = 'circ';         % Microlens aperture shape ['rect' or 'circ']
ulens_fill_factor = 100;        % Microlens array fill factor (percent)
