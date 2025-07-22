 function [Position,Velocity, anom] = orb2pltdta(DIM,SAT,mu)
% This function takes the information of the orbit and outputs the position
% vectors and velocity vectors throughout the orbit for al values of true
% anomally

% Inputs:
%   DIM        - Dimension (2 or 3)
%   SAT        - Keplerian orbital elements 
%   mu         - Standard gravitational parameter of central body [m^3/s^2]
% Output:
%   Position   - Array of positions along orbit [m]
%   Velocity   - Array of velocities along orbit [m/s]
%   anom       - Corresponding true anomaly for the position [deg]
  
%% Input: orb = [a e i Omega w; vo]
  toll = 100;
    a = SAT.a;
    p = SAT.p;
    e = SAT.e;
    i = SAT.i;
    O = SAT.Omega;              % longitude of ascending node
    w = SAT.w;                  % Argument of periapsis
    Vo = SAT.Vo;                % True anomaly
    
%% 3D

    R_Trans = TransMat(i,O,w);
    
%   V_PQW = zeros(360*toll,3);
%   V_IJK = zeros(360*toll,3);
    % R_PQWa = zeros(360*toll,3);
    % R_PQWb = zeros(360*toll,3);
    % R_IJK = zeros(360*toll,3);
    j = 1;

    for k = linspace(1,361,36100) % [deg] values of the True anomaly around orbit
        anom(j) = k;

        % 2D
            % Find the position vector in perifocal 
            %R_PQWa:
            Rvalpa = a*(cosd(anom(j))-e); % [m] orbit positions
            Rvalqa = a*sqrt(1 - e^2)*sind(anom(j)); % [m] orbit positions
              R_PQWa(j,1) = Rvalpa.';
              R_PQWa(j,2) = Rvalqa.';
              R_PQWa(j,3) = 0;

            %R_PQWb:
            rval = p/(1+e*cosd(anom(j)/toll));
            Rvalp = rval*cosd(anom(j)/toll);
            Rvalq = rval*sind(anom(j)/toll);
            R_PQWb(j,1) = Rvalp;
            R_PQWb(j,2) = Rvalq;
            R_PQWb(j,3) = 0;

            % Calculate the velocity vector in the PQW frame
            vval = sqrt(mu / p);
            Vvalp = vval * -sind(anom(j));
            Vvalq = vval * (e + cosd(anom(j)));
            V_PQW(j, 1) = Vvalp;
            V_PQW(j, 2) = Vvalq;
            V_PQW(j, 3) = 0;

        % 3D
            % Find position 
            %R_IJK = R_Trans*R_PQW;
            R_ijk = R_Trans.*(R_PQWa(j,1:3));
            R_IJK(j,1) = R_ijk(1,1)+R_ijk(1,2)+R_ijk(1,3);
            R_IJK(j,2) = R_ijk(2,1)+R_ijk(2,2)+R_ijk(2,3);
            R_IJK(j,3) = R_ijk(3,1)+R_ijk(3,2)+R_ijk(3,3);     

            % Transform the velocity vector to the IJK frame
            V_ijk = R_Trans .* V_PQW(j, 1:3);
            Velocity_IJK(j, 1) = V_ijk(1, 1) + V_ijk(1, 2) + V_ijk(1, 3);
            Velocity_IJK(j, 2) = V_ijk(2, 1) + V_ijk(2, 2) + V_ijk(2, 3);
            Velocity_IJK(j, 3) = V_ijk(3, 1) + V_ijk(3, 2) + V_ijk(3, 3);
            

            j = j+1;

    end
    

             
switch DIM
    case 2
        Position = R_PQWb;
        Velocity = V_PQW;
    case 3
        Position = R_IJK;
        Velocity = Velocity_IJK;
end

% %% Output:
%     switch DIM
%         case 2
% 
%             Orbit = R_PQWb; 
% 
%             %V_PQW:
%               vval1 = sqrt(mu/p);
%                 Vvalp1 = vval1*-sin(Vo);
%                 Vvalq1 = vval1*(e + cos(Vo));
%                   V_PQW1(1) = Vvalp1;
%                   V_PQW1(2) = Vvalq1;
%                   V_PQW1(3) = 0;
%             Velocity = V_PQW1; 
% 
%            %Position:              
%                 Rvalp1 = SAT.r*cos(Vo);
%                 Rvalq1 = SAT.r*sin(Vo);
%                   R1(1) = Rvalp1;
%                   R1(2) = Rvalq1;
%                   R1(3) = 0;
%             Position = R;
% 
%         case 3
%             Orbit = R_IJK;
%             Velocity = SAT.V;  
%             Position = SAT.R;  
%     end
% 


end