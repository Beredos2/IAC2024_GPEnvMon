% --SYM_1--
%Sym_1 is the main executable file for IAC2023 simulation. The simulation
%consists of a femtosatellite swarm in low earth orbit, which is collecting
%data from the world magnetic model. A perturbation is applied to the world
%magentic model, which affects the magnetic field in a scale of minutes,
%seconds, and microseconds. At the end of the simulation, the data is
%processed to return intensity maps over time of a single cell, intensity
%maps over space for multiple cells, and a field reconstruction of a cell
%over time. 

%% Initialize
%Clean up 
clear 
clc
addpath(genpath("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\IAC2024_GPxFemtosatsProject"))
%% Time parameter definition-------------------------------------------------
n_orbits = 100; %Choose the number of orbits we want to analyse
t_0 = 1672531200; %mission start epoch [s] %NOTE: 00:00 1 of January, 2023
t_f = t_0+n_orbits*5.43010002e3; %mission end epoch [s] %NOTE: 00:00 2 OF January, 2023 ; 5.43010002e3 is the orbital period
dt = 1; %stepsize [s]
%Initial orbit design------------------------------------------------------
a = 7400; %Semi major axis [km]
e = 0.0000001; %Eccentricity
i = 0; %Inclination [rads]
Omega = 0; %Right ascention of the ascending node [rad]
w = 0; %Angle of perigee [rad]
M = 0; %Mean anomaly [rad]
Kepler = [a,e,i,Omega,w,M]; %Keplerian vector 
%Swarm definition----------------------------------------------------------
ID = 1; %identification digit
type = "MagMeasure"; %Swarm Type
%Drone manufacturing instructions
dType = "Drone_m";
dPop = 163; 
dNeed = 50; 
systems = table(dType,dPop,dNeed);
%Planetary Body parameters-------------------------------------------------
PRadius = 6371; %[km] planetary radius
rotRate = [0 0 2*pi/86400]; %[rad/s] Rotation rate relative to ECI
R_p = [0 0 0]; 
V_p = [0 0 0];
Theta_p = [0 0 0]; %NOTE: must calculate for the starting epoch
% Instantiation 
myClock = tiktok(t_0,t_f,dt); 
mySwarm = swarm(ID,type,systems,Kepler,t_0);
myEarth = celeBod(R_p,V_p,Theta_p,rotRate,PRadius,myClock); 
%Clean up
clear a dNeed dPop dt dType e i ID Kepler M Omega PRadius R_p rotRate systems t_0 t_f Theta_p type V_p w
disp("instantiation complete")
%__________________________________________________________________________
%% Orbit propagation
deployBrood(mySwarm); %Apply deployment delta vee
for i=1:length(mySwarm.network{1})
    propagateOrb(mySwarm.network{1}{i},myClock);
end
disp("orbit propagation complete")
%__________________________________________________________________________
% %% Data collection 
% %No error, yes perturbation
% tic
% for i=1:length(mySwarm.network{1})
%     getMdata(mySwarm.network{1}{i},myEarth,myClock)
% end
% disp("data collection complete")
% toc
%__________________________________________________________________________ 
runtime = toc;
%
load handel
sound(y,Fs)

clear y Fs 


