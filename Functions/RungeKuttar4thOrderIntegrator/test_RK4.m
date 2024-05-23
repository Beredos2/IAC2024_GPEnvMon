%==========================================================================
%  Satellite Propagation with Runge-Kutta 4th Order Numerical Integration
%==========================================================================
%> @author Ozan Kilic
%> Middle East Technical University, Space Geodesy Division
%> contact: kilic.ozan@metu.edu.tr
%==========================================================================
%> Inputs: Initial Satellite PV Vector in ECI Frame
%
%        Position and Velocity in (x,y,z,vx,vy,vz) (in meters)
% 
%        OR
% 
%        Kepler Orbital Elements:
%        a            = Semi-major axis (km)
%        e            = Eccentricity (unitless)
%        i            = Inclination (rad)
%        Omega        = Right-Ascension of the ascending node (rad)
%        w            = Argument of periapse (rad)
%        M            = Mean anomaly (rad)
%        mu           = Gravitational parameter of the system (GM), km^3/s^2
%
%        h            = Step Size
%        steps        = Number of Steps
% 
%> Output: Propagated Satellite PV Vector in ECI Frame
%==========================================================================
clc
clear
close all
format long g
%==========================================================================
% Initial Coordinates
%==========================================================================
% Cartesian Position and Velocity
% X = [1892.775;  28831.100; 6415.927; -2.140; -0.513;  2.937]*1e3;

% Keplerian Elements
a = 29599.8;
e = 0.0001;
i = 0;
Omega = 0;
w = 0;
M = 0;
[RECI, VECI] = Kepler2RV(a, e, i, Omega, w, M);
X = [RECI;VECI]*1e3;
%==========================================================================
% Propagation
%==========================================================================
h = 300;                  % Step Size (sec)
steps = 300;              % Number of Steps

[X_RK] = RK_4(X,h,steps); % Orbit in ECI (PV)

%=========================================================================