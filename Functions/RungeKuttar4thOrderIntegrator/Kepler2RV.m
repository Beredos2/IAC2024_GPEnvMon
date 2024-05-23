%==========================================================================
%> Kepler2RV
%==========================================================================
%> This function converts keplerian elements to position and velocity
%> vector.
%>
%> @author Ozan Kilic
%> Middle East Technical University, Space Geodesy Division
%> contact: kilic.ozan@metu.edu.tr
%==========================================================================
% Inputs:    a       = Semi-major axis (km)
%            e       = Eccentricity
%            i       = Inclination (rad)
%            Omega   = Right-Ascension of the ascending node (rad)
%            w       = Argument of periapse (rad)
%            M       = Mean Anomaly (rad)
%            mu      = Gravitational parameter of the system (GM)(km^3/s^2)
%
% Outputs:   RECI    = (x,y,z)    Cartesian Position (km)
%            VECI    = (vx,vy,vz) Cartesian Velocity (km/s)
%==========================================================================
function [RECI, VECI] = Kepler2RV(a, e, i, Omega, w, M)

mu  = 398600.441799999971;                 % Earth's Standard Gravitational Parameter (GM)
Eps = 2.22044604925031e-8;                 % Machine epsilon

E   = SolveKepler(M,e,Eps);                % Eccentric Anomaly from Mean Anomaly
M   = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));  % True Anomaly obtained from Eccentric Anomaly
p   = a*(1 - e^2);                         % Semiparameter (km)

%  Position Coordinates in Perifocal Coordinate System

x   = (p*cos(M)) / (1 + e*cos(M));         % x-coordinate (km)
y   = (p*sin(M)) / (1 + e*cos(M));         % y-coordinate (km)
z   = 0;                                   % z-coordinate (km)
vx  = -(mu/p)^(1/2) * sin(M);              % velocity in x (km/s)
vy  = (mu/p)^(1/2) * (e + cos(M));         % velocity in y (km/s)
vz  = 0;                                   % velocity in z (km/s)

%  Transformation Matrix (3 Rotations)

t_rot = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i) ...
        (-1)*cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i) ...
        sin(Omega)*sin(i); ...
        sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i) ...
        (-1)*sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i) ...
        (-1)*cos(Omega)*sin(i); ...
        sin(w)*sin(i)  cos(w)*sin(i)  cos(i)];

%  Transformation Perifocal xyz to ECI

RECI = t_rot*[x y z]';
VECI = t_rot*[vx vy vz]';

end

function E = SolveKepler(M,e,Eps)
% Iterative Solution of Kepler's equation (M = E-e*sin(E))
%==========================================================================
% Inputs:    M   =  Mean Anomaly (rad)
%            e   =  Eccentricity
%            Eps =  Machine epsilon
% Output:    E   =  Eccentric Anomaly (rad)
%==========================================================================
E0  = M;
E1 = E0 - (E0 - e*sin(E0) - M)/(1  -e*cos(E0));
while (abs(E1-E0) > Eps)
    E0 = E1;
    E1 = E0 - (E0 - e*sin(E0) - M)/(1 - e*cos(E0));
end
E = E1;

end