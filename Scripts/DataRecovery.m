%Data recovery: This script recovers data from the MDIPS demosntration. The
% user determines the time step at which we want to sample, and then the
% script recovers data at that rate. The data is organized as a cell array
% with a number of elements given by the total time steps divided by the
% sample rate. Each cell in the array is a time frame, represented by an
% array. Each line of the array contains the location r, magnetic vector M,
% and magnetic intensity norm(M) as [x y z Mx My Mz norm(M)] 

%0. Load simulation data
addpath(genpath("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\IAC2024_GPxFemtosatsProject"))
load("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\IAC2024_GPxFemtosatsProject\Results\Experiment_3") 
%clean up unnecessary variables
clear D D_p deltaMatrix F index mySTS mySTS2 n R r_o runtime t T

%% Create an analyisis timeline 
% Known sample rate of the simulation: 0.5 seconds 
% If I want to sample 1/second, then I need to sample from myClock.timeline
% once every 2 measurements. 
modelDT = 0.5; %seconds
deltaT = 10; %seconds
T_index = deltaT/modelDT; %index rate 
% Timeline of interest 
timeline = myClock.timeline(:,1); 
timeline = timeline(1:T_index:length(timeline),1);
% Create time-wise data structure
Data = cell(1,length(timeline)); 
% Buffer array
data = zeros(length(mySwarm.network{1}),7); 
% Initiate simulated data index counter
T = 1; 
% Assignment loop 
for t = 1:T_index:length(myClock.timeline)
    for i = 1:length(mySwarm.network{1})
        data(i,1:3) = mySwarm.network{1}{i}.r(:,t)';                          % Assign position
        data(i,4:7) = mySwarm.network{1}{i}.D_m(t,:);                         % Assign magnetic field measurement 
    end
    Data{T} = data;  
    T = T+1;                                                         % Assign data to frame
end

clear ans data myEarth mySwarm t i T myClock modelDT deltaT