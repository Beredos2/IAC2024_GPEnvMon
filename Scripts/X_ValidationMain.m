%Cross validation main is a scritp that loads a data set from an MDIPS simulation. 

% All data is arranged into a 1xk cell array, where k is the number of time frames
% being studied. Each cell contains an nx7 array, so that it holds all
% drones, the position [x y z], magnetic vector [Mx My Mz], and magnetic
% vector norm M. 
clc 
clear
addpath(genpath("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\IAC2024_GPxFemtosatsProject"))
%% get data
DataRecovery
D = Data; 
clear Data dataCyl
%*************************************************************************%
%% Partition data p=0.15 c=4
%Initialize parameters
T = 2000; % time of measurement
percentage = 0.15; % percentage of points to be removed for cross validation
iterations = 100; % number of iterations
kernels = ["exponential";"squaredexponential";"matern32";"ardmatern52";"rationalquadratic"];
% assign storage space
MAE = cell(length(kernels),1);
MAPE = cell(length(kernels),1);
RMSE = cell(length(kernels),1);
myTable = table(kernels,MAE,MAPE,RMSE); 
% primary calc loop - NOTE:could use paralel computing to speed up loop
for k = 1:length(kernels)
    %initialize storage vectors
    MAEk = zeros(1,iterations);
    MAPEk = zeros(1,iterations); 
    RMSEk = zeros(1,iterations); 
    for i=1:iterations
        % random partition
        [trainSet,testSet] = partitionArray(percentage,D,T);
        trainSet = sort(trainSet);
        testSet = sort(testSet);
        %train GP model 
        X = D{T}(trainSet,3); %position
        Y = D{T}(trainSet,7); %magnetic field intensity
        myModel = fitrgp(X,Y,'KernelFunction',kernels(k)); 
        yPred = predict(myModel,D{T}(testSet,3)); 
        %collect and store performance data
        [MAEk(i),MAPEk(i),RMSEk(i)] = gpPerformance(D,testSet,yPred,T); %Store into available memory
    end 
    %Calculate statistics
    meanMAE = mean(MAEk);
    stdMAE = std(MAEk); 
    meanMAPE = mean(MAPEk);
    stdMAPE = std(MAPEk); 
    meanRMSE = mean(RMSEk);
    stdRMSE = std(RMSEk); 
    %Update model performance Table
    myTable.MAE{k} = [meanMAE,stdMAE]; 
    myTable.MAPE{k} = [meanMAPE,stdMAPE]; 
    myTable.RMSE{k} = [meanRMSE,stdRMSE]; 
end

%Display table
disp(myTable); 

clear i j MAEk MAPEk RMSEk meanMAE meanMAPE meanRMSE stdMAE stdMAPE stdRMSE

%% -FUNCTIONS- 
%--------------------------------------------------------------------------
% CARTESIAN2POLAR(D)-------------------------------------------------------
% is a function that changes the coordinate system from cartesian to
% cylindrical coordinates. It takes D and Ti as arguments, where D is the
% data set D{1:k}(1:n), and Ti is the time instant of interest. 
function D_pol = cartesian2polar(D,Ti)
    n = height(D{Ti});
    D_pol = zeros(n,7); 
    for i=1:n
        [theta,rho,z] = cart2pol(D{Ti}(i,1),D{Ti}(i,2),D{Ti}(i,3));  
        D_pol(i,1:3) = [theta,rho,z];
        D_pol(i,4:7) = D{Ti}(i,4:7);        
    end 
end
%--------------------------------------------------------------------------
% PARTITIONARRAY(percentage,clusterSize,Data,time_index)-------------------
% Is a function that takes a Data set: Data{1:k}(1:n), looks at the array
% selected by the time_index, and paritions that array into two sets: the
% tranining set "TrSet" and the test set "TstSet". The function user also
% provides the proportion of satellites "p" that must be assigned to the
% testing set.

function [TrSet, TstSet] = partitionArray(p,D,Ti)
% Establish size of data set "n" 
n = height(D{Ti});
% Check that x2 is an integer, and if not then reassign p so that x2
% approaches its desired value from below. 
x2 = n*p;% size of X2
if x2>floor(x2)
    x2 = floor(x2); 
    p = x2/n;
    disp("WARNING: x2 is not an integer, so p has been recalculated")
    disp(p)       
end
x1 = n-x2;
fprintf('X1 has %i elements \n',x1)
fprintf('X2 has %i elements \n',x2)

% ESTABLISH SYSTEM BOUNDARIES AND PROTECT THEM FROM ELIMINATION
% Boundaries defined in cylindrical cooridnates
Dpol = cartesian2polar(D,Ti); 
% Find the boundary satellites in cylindrical coordinates
myADJ = getADJ(Dpol);
%Theta boundaries
[thetaB1,thetaB2] = find(myADJ{1}==max(max(myADJ{1})));
%Rho boundaries
[rhoB1,rhoB2] = find(myADJ{2}==max(max(myADJ{2})));
%Z boundaries
[zB1,zB2] = find(myADJ{3}==max(max(myADJ{3})));
forbiddenList = [thetaB1,thetaB2,rhoB1,rhoB2,zB1,zB2]; 

%Initialize partition variables
X1 = forbiddenList; 
X2 = [];
%Initialize ADJ
myADJ = getADJ(D,Ti); %WARNING: recycling variable name
myADJ = sqrt(myADJ{1}.^2+myADJ{2}.^2+myADJ{3}.^2); %straight line distance between satellites
cCount = 0;

% PARTITION ASSIGNMENT LOOP
while length(X2)<x2 
    % protect boundaries and already assigned points
    forbiddenList = [X1 X2]; 
    % list satellites that can be picked
    allowedList = setdiff(1:n,forbiddenList); 
    % randomly select a cluster nucleus
    target = allowedList(randi(length(allowedList))); 
    % generate a group size "c" within constraints
    if cCount == 0
        %don't allow c=x2 on first iteration
        c = randi(floor(x2/2)); % FIX: ***Force at least size 2***
    else
        %don't allow c to be greater than number of slots left
        c = randi(x2-length(X2)); %FIX: ***See notes***
    end
    % select c-1 nearest neighbours
    killList = find(ismember(myADJ(target,:),mink(myADJ(target,:),c)));
    % assing nearest neighbours + target to X2 and forbiddenList 
    X2 = [X2 killList];
    forbiddenList = [forbiddenList killList];
    cCount = cCount+1;
end
% Assign non-test points to the training set
X1 = setdiff(1:n,X2);
TrSet = X1; 
TstSet = X2; 
end


function [MAE,MAPE,RMSE] = gpPerformance(D,testSet,y_pred_1,T)
    % Calculate stats
    % Root Mean Square Error (RMSE) for the test set
    RMSE = sqrt((1/length(testSet))*sum((y_pred_1 - D{T}(testSet,7)).^2));
    % Mean Absolute Error
    MAE = (1/length(testSet))*sum(abs(y_pred_1 - D{T}(testSet,7)));
    % Mean Aboslute Percentage Error
    MAPE = (1/length(testSet))*sum(abs((y_pred_1 - D{T}(testSet,7))./D{T}(testSet,7)));
end


