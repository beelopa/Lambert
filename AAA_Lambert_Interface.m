clear,clc,close all
format long;
%% Input:
Pram = {'a' 'e' 'i' 'Omega' 'w' 'Vo'};      % Orbital Element Choice for orb Function
Method = [1];                               % 1)Universal Variable Method 2)Lagrange Method .... 
tolerance = 1e-5;                           % Tolerance
Dim = 3;
EPH = '421';
Center = 'SolarSystem';

% Starting and end Time:
T1 = datetime(2020,7,1,0,0,0);       	%[yyyy/mm/dd] Launch Date  
T2 = datetime(2021,12,20);          	%[yyyy/mm/dd] Arrival Date  

% Bodies:
    % Sun:
    Body_Center = bodysort();
    Body_Center.name = 'Sun';
    Body_Center.mass = (1.989e30);   % [kg] Sun Mass
    Body_Center.radius=6.957e8;      % [m]  Sun Radius @ Equator
    % Earth:
    Body_Start = bodysort();
    Body_Start.name = 'Earth';
    Body_Start.mass = 5.9722e24;     % [kg] Earth Mass
    Body_Start.radius = 6378e3;      % [m]  Earth Radius @ Equator
    % Mars:
    Body_Target = bodysort();
    Body_Target.name = 'Mars';
    Body_Target.mass=0.64171e24;  % [kg] Mars Mass
    Body_Target.radius=3396e3;    % [m]  Mars Radius @ Equator

%% Constants:
G = 6.674083131313131313e-11;       % Gravitational Constant [m^3/(kg/s^2)]
     
% Unit Vector:
HAT = struct('i',[1 0 0],'j',[0 1 0],'k',[0 0 1]); 
    
%% INITIAL Values:
u  = @(M) G*M;                         % �
mu = u(Body_Center.mass);             % [m^3/s^2]
     
% Invariable Plane Transfer Matrix:
ICR2_TRANS1 = ICRF2IVP(T1,11,EPH,Center);
ICR2_TRANS2 = ICRF2IVP(T2,11,EPH,Center);
     
% Body Info:
fprintf('Kepler Elements for Planetary Bodies:[%s %s %s %s %s; %s]\n',...
       Pram{(1)},Pram{(2)},Pram{(3)},Pram{(4)},Pram{(5)},Pram{(6)});  
   
for body = ['C' 'S' 'T']
    % Grab Data:
        switch body
            case 'C'
                Body = Body_Center;
            case 'S'
                Body = Body_Start;
            case 'T'
                Body = Body_Target;
        end
         
    % Get Body Info:
        % Position and velocity in ecliptic plane 
        [Body.P1_ecl,Body.V1_ecl] = planetEphemeris(juliandate(T1),Center,Body.name,EPH,'km'); 
        [Body.P2_ecl,Body.V2_ecl] = planetEphemeris(juliandate(T2),Center,Body.name,EPH,'km'); 
         
        % Position and velocity in invariable plane
        Body.P1= matsolv(Body.P1_ecl,ICR2_TRANS1)*10^3; 
        Body.V1= matsolv(Body.V1_ecl,ICR2_TRANS1)*10^3;
        Body.P2= matsolv(Body.P2_ecl,ICR2_TRANS2)*10^3; 
        Body.V2= matsolv(Body.V2_ecl,ICR2_TRANS2)*10^3;
            
        % Find orbital elements corresponding to positions above
        [Body.E6, Body.Initial] = rv2orb(Body.P1, Body.V1, mu); 
        [Body.E6, Body.Final  ] = rv2orb(Body.P2, Body.V2, mu);
        Body.Avg = AvgKepElm(Body.Initial, Body.Final);
         
        % Find positions along the orbit (plotting pruposes):
        [Body.Path1, Body.P11, Body.V11] = orb2pltdta(Dim,Body.Initial,mu);
        [Body.Path2, Body.P21, Body.V21] = orb2pltdta(Dim,Body.Final,mu);         
        [Body.Path3,~,~] = orb2pltdta(Dim,Body.Avg,mu);
         
        % Print information 
        fprintf('  %s:\n',Body.name);  
        fprintf("    Initial Orbit = [%.5g,%.5f,%.3f\260,%.3f\260,%.3f\260;%.3f\260]\n",...
                    Body.Initial.a,    Body.Initial.e,Body.Initial.i,...
                    Body.Initial.Omega,Body.Initial.w,Body.Initial.Vo);
         
        fprintf("    Final Orbit   = [%.5g,%.5f,%.3f\260,%.3f\260,%.3f\260;%.3f\260]\n",...
             Body.Final.a,    Body.Final.e,Body.Final.i,...
             Body.Final.Omega,Body.Final.w,Body.Final.Vo);
         
    % Save Data:
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
                                                                     
% Radius:
R1 = Body_Start.P1;    % Starting Position
r1 = norm(R1);         % Magnitude of Radius 1 (m)
R2 = Body_Target.P2;   % Ending Position
r2 = norm(R2);         % Magnitude of Radius 2 (m)
      
% Time of Flight:
[~,~,~,toh,tom,tos] = datevec(between(T1,T2,'time'));        % Time of Flight in hours
TOF = 60*(60*toh + tom) + tos;                               % Time of Flight in Seconds
N   = round((TOF/Body_Target.Initial.TP)-.5);                % Number of Orbits
    
clear ICR2_TRANS2 u G EPH tom tos toh body Body
            
%% Solver:

% Select method to solve Lambert's problem with
for Eq = Method
    Satellite = bodysort();
      Satellite.P1 = R1;
      Satellite.P2 = R2;
    % Solving via Universal Variable Method
    if Eq == 1          
       Satellite.name      = 'Satellite1';
       Satellite.Method    = 'Universal Variable';
       Satellite.IttMethod = 'Bisection';
    % Solving via Lagrange Method
    elseif Eq == 2       
       Satellite.name      = 'Satellite2';
       Satellite.Method    = 'Lagrange';
       Satellite.IttMethod = 'Bisection';
    % Solving via Gauss Method 
    elseif Eq == 3      
       Satellite.name      =  'Satellite3';
       Satellite.Method    = 'Gauss';
       Satellite.IttMethod = 'Newton';
    % Solving via Giulio Avanzini's Method 
    elseif Eq == 4      
       Satellite.name      =  'Satellite4';
       Satellite.Method    = 'Giulio Avanzini';
       Satellite.IttMethod = 'Bisection';
    end


% Solving Lambert's Problem
[V1,V2,~,run_time,count] = lambert(mu,R1,R2,TOF,N,Satellite.Method,Satellite.IttMethod,tolerance);
% Finding hyperbolic excess velocities at departure and arrival
V1_excess    = norm((V1 - Body_Start.V1)/1000);
V2_excess    = norm((V2 - Body_Target.V2)/1000);
TimeofFlight = (TOF/3600)/24;
       
% Finding orbital elements of Lambert arc
[~,Satellite.Orbit] = rv2orb(R1,V1,mu);
         
% Plot Information:
[Satellite.Path1,Satellite.P11,Satellite.V11] = orb2pltdta(Dim,Satellite.Orbit,mu);
      
% Print:
      fprintf("%s Method:     %s = %f           [%s %s %s %s %s; %s]\n",...
        Satellite.Method, Satellite.Variable,Satellite.Var,...
        Pram{(1)},Pram{(2)},Pram{(3)},Pram{(4)},Pram{(5)},Pram{(6)});
      fprintf("\t%s: Run Time = %f sec     Iterations = %f\n",...
        Satellite.IttMethod,run_time,count);
      fprintf("\t  orbit = [%.5g,%.5f,%.3f\260,%.3f\260,%.3f\260;%.3f\260]\n",...
        Satellite.Orbit.a,Satellite.Orbit.e,Satellite.Orbit.i,...
        Satellite.Orbit.Omega,Satellite.Orbit.w,Satellite.Orbit.Vo);
      fprintf("\t  V1 = [%f,%f,%f]\n",...
        V1(1),V1(2),V1(3));
      fprintf("\t  V2 = [%f,%f,%f]\n",...
        V2(1),V2(2),V2(3));
      fprintf("\t  V1 Hyperbolic = %f\n",...
        V1_excess);
      fprintf("\t  V2 Hyperbolic = %f\n",...
        V2_excess);
      fprintf("\t  Time of Flight = %f\n",...
        TimeofFlight);
      fprintf("\n");
      
% Save Data:
      switch Eq
          case 1
              Satellite1 = Satellite;
          case 2
              Satellite2 = Satellite;
          case 3
              Satellite3 = Satellite;
          case 4
              Satellite4 = Satellite;
      end
      clear Satellite
end
clear count run_time Eq Var

%% Visual:
    % Criteria:
       Planet_Scale = 1500;
            
    % Plot: A View of the Transfer  
       Graph = figure(1);hold on
       title('3D Planet Plot');
         Graph.WindowStyle = 'docked';
         Graph.Renderer = 'opengl';
         Graph.Color = 'white';
         axis tight;box on;grid on;
         xlabel('x [m]');    ylabel('y [m]');    zlabel('z [m]');
         xlim([-1.5*r2 1.4*r2]); ylim([-1.5*r2 1.5*r2]); zlim([-1.5*r2 1.5*r2]);
          
      % Orbital Elements:   
         plot3([0,2e11],[0,0],   [0,0],   'k');
         plot3([0,0],   [0,2e11],[0,0],   'k');
         plot3([0,0],   [0,0],   [0,2e11],'k');
           text(2e11,0,0,'I');text(0,2e11,0,'J');text(0,0,2e11,'K')
      
      % Plotting Bodies & Orbits:    
        % Plot Center Planet: [Sun]
           S_Size = 15*Body_Center.radius; 
             Body_Center.plot_i = plotbody(3,[0 0 0],S_Size,36,[1.00 1.00 0.00]);%        
%              text(Body_Center.P1(1),Body_Center.P1(2),Body_Center.P1(3),'Sun');
             
        % Plot Orbiting Body: [Earth]
           E_Size = Planet_Scale*Body_Start.radius;%/(Sun.radius);           
             Body_Start.plot_i = orb2plt(3,E_Size,[0.75 1.00 0.75],...    %[0.75 1.00 0.75]'.g',Plot_Scale*Earth.radius/(Sun.radius)
                 Body_Start.Path1,Body_Start.P11,Body_Start.V11); 
             text(Body_Start.P11(1),Body_Start.P11(2),Body_Start.P11(3),'Earth:Initial')
             
             Body_Start.plot_f = orb2plt(3,E_Size,[0.75 1.00 0.75],...    %[0.75 1.00 0.75]'.g',Plot_Scale*Earth.radius/(Sun.radius)
                 Body_Start.Path1,Body_Start.P21,Body_Start.V21); 
%              text(Body_Start.P21(1),Body_Start.P21(2),Body_Start.P21(3),'Earth:Final')
             
        % Plot Orbiting Body: [Mars]
           M_Size = Planet_Scale*Body_Target.radius;%/(Sun.radius);
           
             Body_Target.plot_i = orb2plt(3,M_Size,[0.85,0.33,0.10],...   %'.r' rMars*Plot_Scale, Plot_Scale*rMars/(Sun.radius)
                 Body_Target.Path2,Body_Target.P11,Body_Target.V11); 
%              text(Body_Target.P11(1),Body_Target.P11(2),Body_Target.P11(3),'Mars:Initial')
             
             Body_Target.plot_f = orb2plt(3,M_Size,[0.85,0.33,0.10],...   %'.r' rMars*Plot_Scale, Plot_Scale*rMars/(Sun.radius)
                 Body_Target.Path2,Body_Target.P21,Body_Target.V21); 
             text(Body_Target.P21(1),Body_Target.P21(2),Body_Target.P21(3),'Mars:Final')
        
     % Ploting Transfers: 

                Satellite.plot = plot3(Satellite.PathPosition(:,1),...
                                        Satellite.PathPosition(:,2),...
                                        Satellite.PathPosition(:,3),'--');
                Satellite.plot.Color = [1 0 0];
         
        
   plot3([0,R1(1)],[0,R1(2)],[0,R1(3)],'k');
     text((0+R1(1))/2,(0+R1(2))/2, 'R1')  
   plot3([0,R2(1)],[0,R2(2)],[0,R2(3)],'k');
     text((0+R2(1))/2,(0+R2(2))/2, 'R2')    
   plot3([R1(1),R2(1)],[R1(2),R2(2)],[R1(3),R2(3)],'k');
     text((R1(1)+R2(1))/2,(R1(2)+R2(2))/2, 'C')
     view(3)
     
   hold off
 
clear Eq E_Size S_Size M_Size Planet_Scale
         
         
%% Conclusion:
% fprintf("Even though they used the same number of itterations, Lagrange \n");
% fprintf("method took nearly twice as long as Universal VaR1able Method.\n");
% fprintf("This could be caused by the differences in formulas, as Lagrange\n");
% fprintf("put everything into one equation while Universal is seperated into\n");
% fprintf("multiple functions.\n");
%%
%     if xb <= xa || (Fa * Fb) > 0
%       fprintf("Poor Starting Parameters %f",Fa*Fb);
%       return
%     end


%      while abs(TOF - Fx) > t %Fx>(0+t) || Fx<(0-t)
%          if Fx>TOF
%              x = x + n;
%          elseif Fx<TOF
%              x = x - n;
%          else
%              x = x;
%          end
%          Fx = f(x);
%          i = i+1;
%      end





%       y = fy(x);
%        z = fz(x,y);
%         S1 = fS1(x,y,z);
%          Q = fQ(x,S1);
%           TOFg = fTOF(x,z,Q);


%          for body = {Starting_Body,Destination_Body}
%              if strcmp(Starting_Body,body)
%                  s1 = Earth_i;s2 = Earth_f;
%              else
%                  s1 = Mars_i;s2 = Mars_f;
%              end
%              for var = {'H' 'h' 'N' 'n' 'a' 'p' 'EV' 'e' 'i' 'w' 'II' 'T' 'ra' 'rp' 'TP' 'MM' 'Omega'}
%                  body.var = (s1.var + s2.var)/2;
%              end
%          end
% %              for var = {E M Lo uo Vo}
% %                  body.var = s1.var
                 
