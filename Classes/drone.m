 %------------DRONE---------------------------------------------------------
%The drone is an femtosatellite object, which inherits from handle. 
%------------------------Prorpertites--------------------------------------
        %DEFINITION********************************************************
        %ID      %integer that uniquely identifies each drone
        %type    %string that describes drone funciton
        %ASTRODYNAMICS*****************************************************
        % a       %semi-major axis (km)
        % e       %eccentricity
        % i       %inclination (rad)
        % Omega   %right ascention of the ascending node (rad)
        % w       %argument of periapsis (rad)
        % M       %mean anomaly
        %STATE*************************************************************
        % r       %position (m)
        % v       %velocity (m/s)
        %SCIENCE***********************************************************
        % D       %data storing array as [X, Y, Z, H, D, I, F]

%-------------------------Methods------------------------------------------
%CONSTRUCTOR: drone(kepler)
%kepler: the keplerian orbital parameters as [a;e;i;Omega;w;M]

%PROPAGATOR: propagate(D_obj,r_0,v_0,dT)
%Is a Runge-Kuta numerical integrator that propagates the trajectory of the
%spacecraft from its current position (r_0) to the next. It calculates both
%the position and the orbital velocity over the time interval dT. 
%--D_obj: the spacecraft being propagated
%--r_0: position initial condition
%--v_0: velocity initial condition
%--dT: interval of integration

% UPDATE NOTES:
% - Currently the drones don't self update their orbital position
% (i.e., they keep track of r and v, but their keplerian parameters are not
% updated. Fix in update

classdef drone<handle
    properties 
        %-DEFINITION
        ID      %integer that uniquely identifies each drone
        type    %string that describes drone funciton
        %-ASTRODYNAMICS-
        %Keplerian parameters
        a       %semi-major axis (km)
        e       %eccentricity
        i       %inclination (rad)
        Omega   %right ascention of the ascending node (rad)
        w       %argument of periapsis (rad)
        M       %mean anomaly
        %-STATE-
        r       %position (m) throughout time (Orbit propagation occurs in km, but is converted to m)
        v       %velocity (m/s) velocity throughout time  (Orbit propagation occurs in km, but is converted to m)
        %-SCIENCE-
        D_m     %data storing array, collects magnetic field intensity (nT?) 
    end
    
    methods
        %---------------------------CONSTRUCTOR----------------------------
        function D_obj = drone(ID,type,Kepler)
            %Set Keplerian parameters
            D_obj.ID = ID;
            D_obj.type = type; 
            D_obj.a = Kepler(1);
            D_obj.e = Kepler(2); 
            D_obj.i = Kepler(3); 
            D_obj.Omega = Kepler(4);
            D_obj.w = Kepler(5);
            D_obj.M = Kepler(6);
            %Convert Kepler2RV
            [D_obj.r,D_obj.v] = Kepler2RV(Kepler(1),Kepler(2),Kepler(3),Kepler(4),Kepler(5),Kepler(6)); 
        end
           
        %--------------------------ASTRODYANAMICS--------------------------
        %Propagator: propagates the orbit of the satellite throughout the
        %entire simulation, and stores r and v values for the drone's
        %lifecycle in the r and v variables. 
        function propagateOrb(D_obj,C_obj)  
            %prep IC for RK function
            X = [D_obj.r;D_obj.v]*1e3; %NOTE: r_0 and v_0 are input in kilometers and kilometer/s at the moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            %Propagate orbit using RK_4
            [X_RK] = RK_4(X,C_obj.dT,C_obj.nT);
            %Separate position and velocity into object properties
            D_obj.r = X_RK(1:3,:);
            D_obj.v = X_RK(4:6,:);
        end
        %getMdata: function that generates data about the state of the 
        %magnetic field around the earth. 
        function getMdata(D_obj,P_obj,C_obj)%D_obj==Drone;P_obj==Planet;C_obj==symClock
           D_obj.D_m = zeros(length(C_obj.timeline),4);
           for j=1:length(D_obj.D_m)
               %Construct rotation matrix
               R_z = [cos(P_obj.Theta_t(j)) sin(P_obj.Theta_t(j)) 0 ; -sin(P_obj.Theta_t(j)) cos(P_obj.Theta_t(j)) 0; 0 0 1];
               %Get drone position in planetary cartesian
               R_ECEFc = R_z*D_obj.r(:,j);
               %Get drone position in planetary polar coordinates
               [azimuth,elevation,r_abs] = cart2sph(R_ECEFc(1),R_ECEFc(2),R_ECEFc(3)); %Polar coordinates of drone position
               %convert r to height
               h = r_abs-P_obj.radius*10^3;
               R_ECEFp = [azimuth,elevation,h]; %Drone position in planetary coordinates (longitude,latitude,altitude)
               %Get data
               [XYZ,F] = planetMagField(P_obj,R_ECEFp(3),R_ECEFp(2),R_ECEFp(1),C_obj.timeline(j,2));
               D_obj.D_m(j,:) = horzcat(XYZ(1),XYZ(2),XYZ(3),F); 
           end
        end
    end
end


