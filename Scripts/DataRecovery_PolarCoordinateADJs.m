%This is the first of draft scripts writtent for processing the data
%recovered from and MDIPS architecture, for dicussion and use with math
%colleagues. 

%UPDATE NOTES: Currently we are working with the full data set! Can we
%subtract unwanted data to get a lower resolution picture? 

%0. Load simulation data
addpath("\MATLAB\MDIPS Demo\0DegreeInc")
load("\MATLAB\Results\myData_1.mat") 
%clean up unnecessary variables
clear D D_p deltaMatrix F index mySTS mySTS2 n R r_o runtime t T

%% 1. Recover all data, and arrange into the correct format. 
% For the purpose of storing data in simple to access arrays, each time
% frame is stored as a cell in the Data cell array. Then, each Data cell
% contais a nx7 array, which gives [x,y,z,M_x,M_y,M_z,norm(m)]; 
% Create a data storing framework
Data = cell(1,length(myClock.timeline));                                   % high-level, frame by frame storage of all measurements
d = zeros(length(mySwarm.network{1}),7);                                   % data buffer for holding frame data during recovery cycle 
% Time loop 
for t = 1:length(myClock.timeline)
    %id loop 
    for i = 1:length(mySwarm.network{1})
        d(i,1:3) = mySwarm.network{1}{i}.r(:,t)';                          % Assign position
        d(i,4:7) = mySwarm.network{1}{i}.D_m(t,:);                         % Assign magnetic field measurement 
        
    end
    Data{t} = d;                                                           % Assign data to frame
end

clear ans d myEarth mySwarm t

%% 2. Convert position data to spherical coordinates, and store.
SphericalCoordinates = cell(1,length(Data)); 
RinS = zeros(height(Data{1}),3);                                           % Create a buffer for temporary storage
%time loop 
for t = 1:length(Data)
   %id loop
    for id = 1:height(Data{t})
        [RinS(id,1),RinS(id,2),RinS(id,3)] = cart2sph(Data{t}(id,1), ...   % Stored as [azimuth, elevation, r] 
            Data{t}(id,2),Data{t}(id,3)); 
    end 
    SphericalCoordinates{t} = RinS; 
end
%clean up 
clear RinS t id
%% -- All mission data at a chosen temporal resolution -- 
% Because the path of the satellites is approximately circular, the best
% fit to understanding swarm shape as a function of time is to observe
% their relative coordinates in the spherical frame. 

% Establish a timeline and resolution 
SampleTime = myClock.timeline(1:10:length(myClock.timeline),:); 

% 3.1 Produce adjacency matrix using the spherical coordinate vectors
%Create cell arrays for long term storage
ADJ_azimuth = cell(1,length(SampleTime));
ADJ_elevation = cell(1,length(SampleTime));
ADJ_r = cell(1,length(SampleTime));
% Create buffer arrays
adj_az = zeros(height(Data{1}));
adj_el = zeros(height(Data{1}));
adj_r = zeros(height(Data{1}));

% Adjacency calculation cycle 
for t=1:length(SampleTime)
    for i=1:height(Data{1}) % 1:n
        for j=i:height(Data{1}) %i:n
            % **AZIMUTH**
            %Fill top right corner
            adj_az(i,j) = SphericalCoordinates{t}(j,1) - SphericalCoordinates{t}(i,1); % Azimuth ADJ
            %Fill lower left corner
            adj_az(j,i) = -adj_az(i,j);
            % ELEVATION
            %Fill top right corner
            adj_el(i,j) = SphericalCoordinates{t}(j,2) - SphericalCoordinates{t}(i,2); % Elevation ADJ
            %Fill lower left corner
            adj_el(j,i) = -adj_el(i,j);
            % ALTITUDE
            %Fill top right corner
            adj_r(i,j) = SphericalCoordinates{t}(j,3) - SphericalCoordinates{t}(i,3); % Radius ADJ
            %Fill lower left corner
            adj_r(j,i) = -adj_r(i,j);
        end
    end
    % Assign to time frame
    ADJ_azimuth{t} = adj_az;
    ADJ_elevation{t} = adj_el;
    ADJ_r{t} = adj_r;
end

%cleanup 
clear adj_az adj_r adj_el t i j 











