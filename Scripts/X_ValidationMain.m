%Cross validation main is a scritp that loads a data set from an MDIPS simulation. 

% All data is arranged into a 1xk cell array, where k is the number of time frames
% being studied. Each cell contains an nx7 array, so that it holds all
% drones, the position [x y z], magnetic vector [Mx My Mz], and magnetic
% vector norm M. 
addpath(genpath("C:\Users\risto\OneDrive\Documents\PhD\MATLAB\IAC2024_GPxFemtosatsProject"))
%% get data
%load('Results\3n100dt1_DataRetrieved.mat')
%load('Results\Experiment_1.mat')
%*************************************************************************%
load("Experiment_3")
for t=1:length(D)
    D{t} = D{t}';
end
clear D_p deltaMatrix F index my myClock myEarth mySTS mySTS2 mySwarm n R r_o runtime t T
%*************************************************************************%
%% Express location in cylindrical coordinates
Dpol = cartesian2polar(D); 
%% Partition data
T = 5000;
percentage = 0.3; 
clusterSize = 3;
[trainSet,testSet] = partitionArray(percentage, clusterSize, D, Dpol,T);
clear percentage clusterSize
%% Plot the two sets in different colours
figure(1)
hold on
title("Test vs train sest")
for i=1:length(trainSet)
    j = trainSet(i);
    scatter3(D{2500}(j,1),D{2500}(j,2),D{T}(j,3),'r.')
end

for i=1:length(testSet)
    j = testSet(i);
    scatter3(D{2500}(j,1),D{2500}(j,2),D{T}(j,3),'g.')
end


%% -FUNCTIONS- 
%--------------------------------------------------------------------------
% CARTESIAN2POLAR(D)-------------------------------------------------------
% is a function that changes the coordinate system from cartesian to
% cylindrical coordinates. As an argument, it takes the data structure
% provided as D{1:k}(1:n). 
function D_pol = cartesian2polar(D)
    k = length(D);
    n = height(D{1});
    D_pol = cell(1,k);
    for t=1:k
        d_pol = zeros(n,3)%7); 
        for i=1:n
            [theta,rho,z] = cart2pol(D{t}(i,1),D{t}(i,2),D{t}(i,3));  
            d_pol(i,1:3) = [theta,rho,z];
            %d_pol(i,4:7) = D{t}(i,4:7);        
        end
        D_pol{t} = d_pol; 
    end

end
%--------------------------------------------------------------------------
% PARTITIONARRAY(percentage,clusterSize,Data,time_index)-------------------
% Is a function that takes a Data set: Data{1:k}(1:n), looks at the array
% selected by the time_index, and paritions that array into two sets: the
% tranining set "TrSet" and the test set "TstSet". 

%[Describe how the partition algorithm works]

function [TrSet, TstSet] = partitionArray(percentage,clusterSize,D,Dpol,Ti)

% CHECK AND/OR RECALCULATE PERCENTAGE AND CLUSTER SIZE TO FIT THE SET SIZE
% Get the sizes [x1,x2] of the array partition
x2 = height(D{Ti})*percentage;
portionFlag = 0;
% Apply bias to force integer set size
if x2>floor(x2)
    x2 = floor(x2); 
    portionFlag = 1;
end
% Calculate x1 from new x2 
x1 = height(D{Ti}) - x2;
% Check flag and provide warning
if portionFlag == 1
    % Recalcualte percentage
    p = x2/(x1+x2);
    disp("percentage changed to")
    disp(p);
    disp("to make partition size integer")
else 
    p = percentage;
end
% Get the number of groups [m] from x2 and cluster size
m = x2/clusterSize; 
clusterFlag = 0;
% Apply bias to force integer
if m>floor(m)
    m = floor(m);
    clusterFlag = 1;
end
%Check flag, warn, and provde residue group
if clusterFlag == 1
    m_res = x2 - clusterSize*m;
    disp("last group is not big enough to account for a full cluster, so last group size reduced to")
    disp(m_res)
else
    m_res = 0;
end
%Store group sizes 
M = [m,m_res]; 


% ESTABLISH SYSTEM BOUNDARIES AND PROTECT THEM FROM ELIMINATION
% Find the boundary satellites in cylindrical coordinates
myADJ = getADJ(Dpol,Ti);
%Theta boundaries
[thetaB1,thetaB2] = find(myADJ{1}==max(max(myADJ{1})));
%Rho boundaries
[rhoB1,rhoB2] = find(myADJ{2}==max(max(myADJ{2})));
%Z boundaries
[zB1,zB2] = find(myADJ{3}==max(max(myADJ{3})));
safePoints = [thetaB1,thetaB2,rhoB1,rhoB2,zB1,zB2]; 

%Initialize partitions
X1 = safePoints; 
X2 = [];
%Initialize ADJ
myADJ = getADJ(D,Ti); %WARNING: recycling variable name
myADJ = sqrt(myADJ{1}.^2+myADJ{2}.^2+myADJ{3}.^2); %straight line distance between satellites


% ASSIGN TEST SET (KILL LOOP)
for i=1:floor(x2/M(1))
    forbiddenList = [X1 X2];
    % Remove safe points and already used points
    allowedList = setdiff(1:(x1+x2),forbiddenList);
    %-randomly identify a target 
    target = allowedList(randi(length(allowedList)));
    %-identify a number of nearby neighbours equal to the cluster size
    killList = find(ismember(myADJ(target,:),mink(myADJ(target,:),M(1))));
    %-add set to X2 
    X2 = [X2 killList];
end
% Catch any residual satellites 
if clusterFlag == 1
    forbiddenList = [X1 X2];
    % Remove safe points and already used points
    allowedList = setdiff(1:(x1+x2),forbiddenList);
    %-randomly identify a target 
    target = allowedList(randi(length(allowedList)));
    %-identify a number of nearby neighbours equal to the cluster size
    killList = find(ismember(myADJ(target,:),mink(myADJ(target,:),M(2))));
    %-add set to X2 
    X2 = [X2 killList];
end
% Assign non-test points to the training set
X1 = setdiff(1:(x1+x2),X2);

TrSet = X1; 
TstSet = X2; 
end

