clear,clc,close all
% Inputs 
mu                  = 1.327e20;
T1                  = datetime(2017,7,1,0,0,0);
T2                  = datetime(2019,12,20);
tolerance           = 1e-5;
Satellite           = bodysort();
Satellite.name      = 'Satellite1';
Satellite.Method    = 'Universal Variable';
Satellite.IttMethod = 'Bisection';
Body_Center         = bodysort();
Body_Center.name    = 'Sun';
Body_Center.mass    = 1.989e30;     % [kg] Sun Mass
Body_Center.radius  = 6.957e8;      % [m]  Sun Radius @ Equator
% Earth:
Body_Start          = bodysort();
Body_Start.name     = 'Earth';
Body_Start.mass     = 5.9722e24;      % [kg] Earth Mass
Body_Start.radius   = 6378e3;         % [m]  Earth Radius @ Equator
% Mars:
Body_Target         = bodysort();
Body_Target.name    = 'Mars';
Body_Target.mass    = 0.64171e24;   % [kg] Mars Mass
Body_Target.radius  = 3396e3;       % [m]  Mars Radius @ Equator

% Obtaining original orbit
[Body_Center, Body_Start, Body_Target, TOF, V1, V2, Satellite, ~] = OriginalOrbit(T1, T2, '421', 'SolarSystem', Body_Center, Body_Target, Body_Start, mu, 3, Satellite, tolerance, 1);


% Finding orbital parameters for propagation 
[~,Orbit_base]  = rv2orb(Body_Start.P1, V1, mu);      % Initial Orbit parameters (to find TA, other orbital elements should remain the same)
[~,Orbit_base2] = rv2orb(Body_Target.P2, V2, mu);     % FInal Orbital parameters (to find TA)

% Use these elements for propagation, not assumed values 
a     = Orbit_base.a;
e     = Orbit_base.e;
i     = Orbit_base.i;
Omega = Orbit_base.Omega;
w     = Orbit_base.w;

% Initialise 
j  = 1;

for k = linspace(Orbit_base.Vo,Orbit_base2.Vo,1000) % [deg] values of the True anomaly around orbit
        anom(j) = k;

        [Position_Propagated(j,:), Velocity_Propagated(j,:)] = orb2rv(a, e, i, Omega, w, anom(j), mu);
            
        disp(j)
        j = j+1;
end


%%
figure
plot3(Satellite.PathPosition(:,1), Satellite.PathPosition(:,2), Satellite.PathPosition(:,3), 'k', 'LineWidth', 1);
hold on;
% Plot the perturbed path with a dashed red line
plot3(Position_Propagated(:,1),Position_Propagated(:,2),Position_Propagated(:,3),'r--', 'LineWidth', 2);

% % Plot the points with more professional colors and sizes
plot3(Body_Start.P1(:,1), Body_Start.P1(:,2), Body_Start.P1(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', [0 0.4470 0.7410]); % Blue
plot3(Body_Target.P2(:,1), Body_Target.P2(:,2), Body_Target.P2(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', [0.8500 0.3250 0.0980]); % Red
plot3(Body_Center.P1(:,1), Body_Center.P1(:,2), Body_Center.P1(:,3), 'o', 'MarkerSize', 16, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor', [0.9290 0.6940 0.1250]); % Yellow
legend({'Original', 'Propagated', 'Earth (@T_1)', 'Mars (@T_2)', 'Sun'}, 'Location', 'best', 'FontSize', 16)










function [Body_Center, Body_Start, Body_Target, TOF, V1, V2, Satellite, anomally] = OriginalOrbit(T1, T2, EPH, Center, Body_Center, Body_Target, Body_Start, mu, Dim, Satellite, tolerance, N)

% Function that outputs the original orbit (i.e. solves Lamberts problem
% and gives the orbital elements of the generated orbit
% Inputs:
%   T1          - Time of departure
%   T2          - Time of arrival
%   EPH 
%   Center      - Center of system 
%   Body_Center - Center of transfer
%   Body_Target - Target Body
%   Body_Start  - Start body
%   mu          - Standard gravitational parameter
%   Dim         - Dimension
%   Satellite
% Outputs:
%   Body_Center - Data of central body 
%   Body_Start  - Data of starting body
%   Body_Target - Data of target body
%   TOF         - Time of flight
%   V1          - Initial Hyperbolic excess velocity
%   V2          - Final hyperbolic excess velocity
%   Satellite   - Structure including all the spacecraft's parameters
%   anomally    - true anomally

% Invariable Plane Transfer Matrix:
ICR2_TRANS1 = ICRF2IVP(T1,11,EPH,Center);
ICR2_TRANS2 = ICRF2IVP(T2,11,EPH,Center);

% Finding the position of planets of interest 

for body = ['C' 'S' 'T']
    %Grab Data:
    switch body
        case 'C'
            Body = Body_Center;
        case 'S'
            Body = Body_Start;
        case 'T'
            Body = Body_Target;
    end
         
    % Find position and velocity of planets at departure and arrival in ecliptic plane:
    [Body.P1_ecl,Body.V1_ecl] = planetEphemeris(juliandate(T1),Center,Body.name,EPH,'km'); 
    [Body.P2_ecl,Body.V2_ecl] = planetEphemeris(juliandate(T2),Center,Body.name,EPH,'km'); 
         
    % Find position and velocities in invariable plane
    Body.P1= matsolv(Body.P1_ecl,ICR2_TRANS1)*10^3; 
    Body.V1= matsolv(Body.V1_ecl,ICR2_TRANS1)*10^3;
    Body.P2= matsolv(Body.P2_ecl,ICR2_TRANS1)*10^3; 
    Body.V2= matsolv(Body.V2_ecl,ICR2_TRANS1)*10^3;
          
    % Find orbital elements at start and finish
    [Body.E6, Body.Initial] = rv2orb(Body.P1, Body.V1, mu); 
    [Body.E6, Body.Final  ] = rv2orb(Body.P2, Body.V2, mu);
    Body.Avg                = AvgKepElm(Body.Initial, Body.Final);
         
     %Plot Elements:
     [Body.Path1, Body.P11, Body.V11] = orb2pltdta(Dim,Body.Initial,mu);
     [Body.Path2,Body.P21,Body.V21]   = orb2pltdta(Dim,Body.Final,mu);
     %[Body.Path3,~,~]                 = orb2pltdta(Dim,Body.Avg,    Body.radius,mu);
         
       %Save Data:
        switch body
           case 'C'
              Body_Center = Body;
           case 'S'
              Body_Start = Body;
           case 'T'
              Body_Target= Body;
        end
         
end
fprintf('\n\n');
                                                                     
%Radius:
R1 = Body_Start.P1;    %Starting Position
R2 = Body_Target.P2;   %Ending Position
      
%Time of Flight:
[~,~,~,toh,tom,tos] = datevec(between(T1,T2,'time'));          % Time of Flight in hours
TOF                 = 60*(60*toh + tom) + tos;                 % Time of Flight in Seconds
N                   = round((TOF/Body_Target.Initial.TP)-.5);  % Number of Orbits    

% Solving Lambert's Problem:
[V1,V2,~,~,~]   = lambert(mu,R1,R2,TOF,N,Satellite.Method,Satellite.IttMethod,tolerance);

% Finding orbit of satellite:
[~,Satellite.Orbit] = rv2orb(R1,V1,mu);        
         
% Position and velocity along the orbit for all the values of true anomally 
[Satellite.PathPosition,Satellite.PathVelocity, anomally] = orb2pltdta(Dim,Satellite.Orbit,mu);

end


% ---------------- function 5 ---------------- %
function [P_SRP, F_SRP,F_SRP_vec, a_SRP, Acc_Vector] = solaracceleration(Position, Radiation)
% Find the position and velocity vector in the perifocal reference frame and invarialble
% reference frame from the keplerian elements at a given true anomaly
% Inputs:
%   S_0                   - Solar Constant at 1AU [W/m^2]
%   c                     - Speed of light [m/s]
%   AU                    - 1 AU [m]
%   Distance_Sun          - Semi-major axis [m]
%   shadow                - Shadow coefficient 
%   C_R                   - Coefficient of reflectivity 
%   A_S                   - Absorbing area [m^s]
%   m                     - Mass [kg]
%   Rad_UV                - Unit vector 

% Output:
%   P_SRP                 - Solar Radiation Pressure
%   F_SRP                 - Force due to Solar Radiation Pressure
%   a_SRP                 - Acceleration caused by the force
%   Acc_Vector            - Acceleration vector acting radially outwards [m/s^2]
    

%Repalce "Radiation" in the input of function with the constants you want
%for this function (Satellite related or soalr radiation pressure related)
    % Calculate the distance from the Sun
    Distance_Sun = sqrt(Position(1)^2 + Position(2)^2 + Position(3)^2);

    % Calculate the radial unit vector
    Rad_UV = [Position(1) / Distance_Sun; Position(2) / Distance_Sun; Position(3) / Distance_Sun]';    

    P_SRP           = (Radiation.S_0/Radiation.c) * ((Radiation.AU/Distance_Sun)^2);    % Solar Radiation Pressure 
    F_SRP           = Radiation.shadow * P_SRP * Radiation.C_R * Radiation.A_S;         % Force due to Solar Radiation Pressure
    F_SRP_vec       = F_SRP * Rad_UV;
    a_SRP           = F_SRP/Radiation.m; 
    Acc_Vector      = Rad_UV * a_SRP;                     % Acceleration caused by F_SRP

end 