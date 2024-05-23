%% MissionData: is a script that will review the MDIPS_Demo results in 
% detail to ensure that the results are relevant to the scientific
% objectives. Easy to read graphs, figures, and tables are presented that
% address the mission history and the data collected throughout the
% mission.

%Recover simulation data
%load("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\MDIPS Demo\0DegreeInc\RESULTS_2.mat")
%Provide access to classes and function toolbox
addpath("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\MDIPS Demo\0DegreeInc")
load('MyLatestData.mat')
% 
%% the path of an arbitrary satellite througout the mission

%t = myClock.timeline(:,1);
figure(1)
hold on
title("Path of an arbitrary satellite throughout the misison")
xlabel("X_ECI")
ylabel("Y_ECI")
zlabel("Z_ECI")
for id=1:length(mySwarm.network{1})
    r_t = mySwarm.network{1}{id}.r; 
    scatter3(r_t(1,:),r_t(2,:),r_t(3,:),10,'red','filled','o')
end
hold off
view([1,0,0])
clear r_t

%% FIGURE 2: SINGLE SPACECRAFT POV
figure(2)
hold on
ylabel("Magnetic field intensity [nT]")
xlabel("Time since deployment [h]")
%xlim([0 ])
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{10}.D_m(:,3),'b')
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{10}.D_m(:,2),'g')
plot(myClock.timeline(:,1)/3600,mySwarm.network{1}{10}.D_m(:,1),'r')
legend('Z_{ECI}','Y_{ECI}','X_{ECI}')
% set(gca,'looseInset',get(gca,'TightInset'))
% saveas(gca,'Results/singleSatFullDataSet.jpeg')
hold off


%% Plotting relative distance between satellites over time
%In this section, the swarm boundaries are plotted over time by identifying
%the greatest distnace between any two satellites at any given time. The
%cartesian coordinates are not adequate for the circular path of the
%satellites, so spherical coordinates are used instead. 
% Azimuth: The angle between the cardinal direction an the position of the
% satellite. 
% Elevation: The angle between the equatorial plane (which includes the
% cardianl) and the satellite. 
% Altitude: The distance between the center of the Earth and the satellite 

% --CARTESIAN COORDINATE RETRIEVAL--
%Initiate variables
T_1 = myClock.timeline(1:10:length(myClock.timeline),1); %Build a timeline with only the steps of interest
D = cell(length(T_1),1); %Cartesian coordiante time frames (k frames going down)
R_i = zeros(3,length(mySwarm.network{1}));  %Position buffer (3 lines, n columns)
%Main data recovery loop 
for t=1:length(T_1) %Time loop
    for id = 1:length(mySwarm.network{1}) %Network loop
        R_i(:,id) = mySwarm.network{1}{id}.r(:,t);
    end
    D{t} = R_i;
end
clear R_i t id ans %cleanup step 

% --POLAR COORDINATE RETRIEVAL--
%Conversion to polar coordinates
%Stored as a vector {azimuth; elevation; r}
% Initiate variables
D_p = cell(length(T_1),1); % Polar coordinate time frames
R_pi = zeros(3,length(mySwarm.network{1})); %Polar position buffer
%Main data transformation loop 
for t=1:length(T_1) %time loop 
    for id = 1:length(mySwarm.network{1}) %Network loop
      [azimuth,elevation,r] = cart2sph(D{t}(1,id),D{t}(2,id),D{t}(3,id));
      R_pi(1,id) = azimuth;
      R_pi(2,id) = elevation;
      R_pi(3,id) = r;
    end
    D_p{t} = R_pi; % Assign data to time frame
end
% cleanup step
clear ans id R_pi t azimuth elevation r

%--POLAR ADJACENCY MATRIX CONSTRUCTION--
% The relationship between elements of the swarm is calculated. The
% matrices constructed are weighed adjacency matrices. The calculation for
% a weight follows this standard: The reference object is called the
% referee, and the object relative to the referee is called the referred
% object. The relationship between the referred and the referree is equal
% to the value of the referred minus the value of the referee. This is to
% make the standard follow the same rules for vector subtraction. 
%1. Initiate variables 
ADJ_p = cell(3,length(T_1)); %ADJ_p stores the adjacency matrix for each time frame and for each variable type:
%ADK_p = {ADJ_az1,ADJaz2,...,ADJazk;      % Azimuth envelope timeline
%         ADJ_el1,ADJ_el2,...,ADJ_elk;    % Elevation envelope timeline
%         ADJ_h1,ADJ_h2,..ADJ_hk}         % Height envelope timeline (measured to the center of the Earth)
ADJ_azbuff = zeros(length(mySwarm.network{1})); %Buffer azimuth adjacency matrix
ADJ_elbuff = zeros(length(mySwarm.network{1})); %Buffer elevation adjacency matrix
ADJ_hbuff = zeros(length(mySwarm.network{1})); %Buffer altitude adjacency matrix
%2. Get and assign data
for t=1:length(T_1)
    for i=1:length(mySwarm.network{1}) % 1:n
        for j=i:length(mySwarm.network{1}) %i:n
            % AZIMUTH
            %Fill top right corner
            ADJ_azbuff(i,j) = D_p{t}(1,j) - D_p{t}(1,i); % Azimuth ADJ
            %Fill lower left corner
            ADJ_azbuff(j,i) = -ADJ_azbuff(i,j);
            % ELEVATION
            %Fill top right corner
            ADJ_elbuff(i,j) = D_p{t}(2,j) - D_p{t}(2,i); % Azimuth ADJ
            %Fill lower left corner
            ADJ_elbuff(j,i) = -ADJ_elbuff(i,j);
            % ALTITUDE
            %Fill top right corner
            ADJ_hbuff(i,j) = D_p{t}(3,j) - D_p{t}(3,i); % Azimuth ADJ
            %Fill lower left corner
            ADJ_hbuff(j,i) = -ADJ_hbuff(i,j);
        end
    end
    % Assign to time frame
    ADJ_p{1,t} = ADJ_azbuff;
    ADJ_p{2,t} = ADJ_elbuff;
    ADJ_p{3,t} = ADJ_hbuff;
end
%cleanup step
clear ADJ_azbuff ADJ_elbuff ADJ_hbuff i j t

%% -- GET MAXIMUM RELATIVE VARIABLE VECTORS -- 
% Initiate variables
delAzimuth = zeros(1,length(T_1)); 
delElevation = zeros(1,length(T_1));
delAltitude = zeros(1,length(T_1)); 
%Assignment loop 
for t=1:length(T_1)
    delAzimuth(t) = max(max(ADJ_p{1,t})); 
    delElevation(t) = max(max(ADJ_p{2,t}));
    delAltitude(t) = max(max(ADJ_p{3,t})); 
end
deltaMatrix = [delAzimuth;delElevation;delAltitude]; 
%cleanup 
clear t delAzimuth delElevation delAltitude

% -- PLOT RELATIVE VARIABLES -- 
figure(10)
hold on
title('Azimuth envelope as a function of time')
xlabel('time [h]')
ylabel('\Delta azimuth [degrees]')
plot(T_1/3600,deltaMatrix(1,:)*180/pi,'.r')
%plot(T/3600,deltaMatrix(2,:)*180/pi,'.g')
%yyaxis right
%ylabel('\Delta r [km]')
%plot(T/3600,deltaMatrix(3,:)/1000,'.b')
%legend('\Delta Azimuth', '\Delta Elevation', '\Delta r'); 
hold off

figure(11)
hold on
title('timewise elevation envelope')
xlabel('time [h]')
ylabel('\Delta elevation[degrees]')
plot(T_1/3600,deltaMatrix(2,:)*180/pi,'.r')
hold off

figure(12)
hold on
title('timewise altitdue envelope')
xlabel('time [h]')
ylabel('\Delta h [km]')
plot(T_1/3600,deltaMatrix(3,:)/1000,'.r')
hold off

%% Boundary check conditions 
% Colin has mentioned the that the way the altitdue and azimuth variables
% increase in an apparently unbounded way is problematic, and may point to
% an error in the model or data collection. He has requested verificaiton
% by plotting in other ways. 

% 1. Plot the altitude of all satellites relative to the original altitude. 
%1.a)Constant altitude for reference 
r_o = norm(mySwarm.network{1}{1}.r(:,1));
%1.b)The set of r's to be compared 
R = zeros(163,length(T_1)); 
for n=1:length(mySwarm.network{1})
    for t=1:length(T_1)
        R(n,t) = norm(mySwarm.network{1}{n}.r(:,t));
    end
end
%1.c)Plot all
figure(13)
hold on
plot(T_1/3600,r_o/1000,'ro');
for n=1:length(mySwarm.network{1})
    plot(T_1/3600,R(n,:)/1000,'b')
end
xlabel('time [h]')
ylabel('a [km]')

%% TFSV analyisis
% The purpose of this anlayisis is to observe locally defined features in
% the magentic field. According to the equation used to modify the magnetic
% field, these variations occur locally relative to the latitude.
% Furthermore, given that the Z comopnent experiences the greatest
% perturbation, any locally determined measurement is going to be more
% clearly differentiated in the z component of the magnetic field. The 7th
% hour is chosen for observation, as it exhibits nearly the greatest
% latitude separtion between drones. 

%1. Choose the time 
T_2 = 18000*2; %first column is used for indexation reasons
%2. Sort the data [location, magnitude]
G_mi = zeros(length(mySwarm.network{1}),4);

for i=1:length(mySwarm.network{1})
    G_mi(i,1) = mySwarm.network{1}{i}.r(1,T_2);
    G_mi(i,2) = mySwarm.network{1}{i}.r(2,T_2);
    G_mi(i,3) = mySwarm.network{1}{i}.r(3,T_2);
    G_mi(i,4) = mySwarm.network{1}{i}.D_m(T_2,4);
end

figure(21)
scatter3(G_mi(:,1),G_mi(:,2),G_mi(:,3),10,G_mi(:,4),'filled')
colorbar
set(gca,'colorscale','log')
xlabel('X_{ECI} [m]')
ylabel('Y_{ECI} [m]')
zlabel('Z_{ECI} [m]')

%% FIGURE 22: Stereotaxic interpretation of Figure 4, for a lattice side of 10km. 
%Get the position at which each measurement was made
R = zeros(3,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    R(:,i) = mySwarm.network{1}{i}.r(:,T_2);
end
%Get the magnetic field intensity for each measurement
F = zeros(1,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    F(1,i) = mySwarm.network{1}{i}.D_m(T_2,4);
end
%Get index location for each measurement
mySTS2 = stereoTaxicSpace(5000);
index = zeros(3,length(mySwarm.network{1}));
for i=1:length(mySwarm.network{1})
    index(:,i) = getIndex(R(:,i),1,1,1,mySTS2);
end
%Create cell array analogous to the final heat map. Each cell in the array
%is a cell in the stereotaxic space, and only cells within the maximum and
%minimum indexed positions are aknowledged (to avoid mapping all of space).
Map = cell(max(index(1,:))-min(index(1,:)),max(index(2,:))-min(index(2,:)));
%Assign all measurements made to a cell corresoponding its indexation
for i=1:length(index)
    if index(1,i)-min(index(1,:)) == 0 && max(index(2,:))-index(2,i) > 0
        Map{index(1,i)-min(index(1,:))+1,max(index(2,:))-index(2,i)} = horzcat(Map{index(1,i)-min(index(1,:))+1,max(index(2,:))-index(2,i)},F(i)); 
    end
    if index(1,i)-min(index(1,:)) > 0 && max(index(2,:))-index(2,i) == 0
        Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1} = horzcat(Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)+1},F(i));
    end
    if index(1,i)-min(index(1,:)) > 0 && max(index(2,:))-index(2,i) > 0
        Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)} = horzcat(Map{index(1,i)-min(index(1,:)),max(index(2,:))-index(2,i)},F(i));
    end %EDIT WARNING: Although previously stable, an error started to occur on the second index of the Map variable. To correct, a +1 was removed. 
end
% Homogenize map and convert to heatmap (i.e., apply an averaging or
% interpolation rule to each cell to resolve each cell into a single value
for i=1:width(Map)
    for j = 1:height(Map)
        if Map{j,i} ~= 0
           Map{j,i} = sum(Map{j,i})/length(Map{j,i});
        else
           Map{j,i} = 0;
        end
    end
end
%Correct orientation
Map = flipud(fliplr(Map).');
%Convert cell array into a matrix array and generate heatmap
Map = cell2mat(Map);
figure(22)
heatmap(Map)
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));


clear ans F G_mi i index j Map Ax n_orbits R T_1 r_s theta_s theta_se 
%% Single cell monitoring
%  Monitor the value of a singel cell over time 
% 1. Choose a grid size and create space
mySTS = stereoTaxicSpace(10000); %The volume chosen reprsents a 100km cube. This might be too big when considering real world applications, but is a perfect excample of continous monitoring when compared against the discontinuity of the propagation modelled here
% 2. Identify a populated cell at 1/4 period
R = mySwarm.network{1}{1}.r(:,10);
index = getIndex(R,1,1,1,mySTS); 
% 3. Monitor cell
[~,F] = monitorCel_F(myClock,mySwarm.network{1},mySTS,index);

disp("monitoring complete")


%% Single voxel value over time / fixed space value over time
figure(31)
hold on 
grid on
T_1 = myClock.timeline(:,1);
plot(T_1,F)
xlabel("time [h]")
ylabel("Magnetic field intensity [nT]")
%xlim([1672531200 1.672558350500100e+09])
hold off


