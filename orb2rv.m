function [R, V] = orb2rv(a, e, i, Omega, w, TA, u)
%% Orbital Elements to Position & Velocity Vectors
% Inputs:
%   a     - Semi-major axis [m]
%   e     - Eccentricity [-]
%   i     - Inclination [deg]
%   Omega - Right Ascension of Ascending Node [deg]
%   w     - Argument of Periapsis [deg]
%   TA    - True Anomaly [deg]
%   u     - Standard gravitational parameter [m^3/s^2]
%
% Outputs:
%   R     - Position vector [m]
%   V     - Velocity vector [m/s]

% Calculate time difference
% dt = calculateTimeBetweenAnomalies(TA1, TA2, a, e, u);
%% Convert angles to radians
i_rad     = deg2rad(i);
Omega_rad = deg2rad(Omega);
w_rad     = deg2rad(w);
TA_rad    = deg2rad(TA);

%% Calculate position and velocity in orbital plane (perifocal frame)
% Semi-latus rectum
p = a * (1 - e^2); % [m]

% Distance from focus
r = p / (1 + e * cos(TA_rad));%[m]

% Position in perifocal frame
r_peri = [r * cos(TA_rad); r * sin(TA_rad); 0];%[m]

% Velocity in perifocal frame
h = sqrt(u * p);  % Specific angular momentum [m]
v_peri = [-sqrt(u/p) * sin(TA_rad); sqrt(u/p) * (e + cos(TA_rad)); 0]; %[m/s]

%% Transformation matrices
Q = TransMat(i_rad, Omega_rad, w_rad);
%% Transform to inertial frame
R = Q * r_peri;
V = Q * v_peri;

% Convert to row vectors (to match your rv2orb format)
R = R';
V = V';

end
