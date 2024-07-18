% Analsis of GP model with time and space as independent variables: this
% script will load mission data via the DataRecovery script, which results
% in a set D{1:k}(1:n,1:7). Then, the script must model the Gaussian
% Process regression for each time frame and carry out a cross validation.
% With each frame evaluated, the error analyisis data is stored into a set
% of arrays that represent MAE, MAPE, and RMSE, as a function of time. The
% final timewise arrays are printed to a plot or set of plots used for
% understanding the method's perfrmance. 

%% initiate variables
% DataRecovery %Check script for dT
% function arguments and loop parameters
p = 0.15; % proportion of removed data points
iterations = 50; %number of iterations used to establish average error anaylisis variables
k = length(Data); % number of time instants %fix applied, previously set to DataLength-1
kernel = "ardsquaredexponential";
% error analysis storage variables: L1 is mean, L2 is std, columns are t
MAE = zeros(2,k-1); 
MAPE = zeros(2,k-1); 
RMSE = zeros(2,k-1);
%% GP_Loop
%Time loop 
for t=2:k 
    % Initiate buffers 
    MAEi = zeros(1,iterations); 
    MAPEi = zeros(1,iterations); 
    RMSEi = zeros(1,iterations); 
    % Iteration loop
    for i=1:iterations
        % Partition array
        [TrSet, TstSet] = partitionArray(p,Data,t);
        TrSet = sort(TrSet); 
        TstSet = sort(TstSet);
        %train GP model 
        X = Data{t}(TrSet,1:3); %position
        Y = Data{t}(TrSet,7); %magnetic field intensity
        myModel = fitrgp(X,Y,'KernelFunction',kernel,'Standardize',1); 
        yPred = predict(myModel,Data{t}(TstSet,1:3)); 
        %collect and store performance data
        [MAEi(i),MAPEi(i),RMSEi(i)] = gpPerformance(Data,TstSet,yPred,t); %Store into available memory
    end
    % Get stats on model performance 
    meanMAE = mean(MAEi);
    stdMAE = std(MAEi); 
    meanMAPE = mean(MAPEi);
    stdMAPE = std(MAPEi); 
    meanRMSE = mean(RMSEi);
    stdRMSE = std(RMSEi); 
    % Assign to final storage
    MAE(1,t) = meanMAE; 
    MAE(2,t) = stdMAE;
    MAPE(1,t) = meanMAPE; 
    MAPE(2,t) = stdMAPE; 
    RMSE(1,t) = meanRMSE; 
    RMSE(2,t) = stdRMSE; 
    disp(t)
end

clear i k kernel iterations p stdMAE stdMAPE stdRMSE MAEi MAPEi RMSEi meanMAE meanMAPE meanRMSE myModel t TrSet TstSet X Y yPred
%% Plot results to screen 

%timeline = timeline0(2:length(timeline0))/3600; 

% figure(1)
% hold on
% plot(MAE(1,:))
% plot(RMSE(1,:))
% 
% legend('MAE','RMSE')

% figure(3)
% title('Correlation betwen MAPE and density along z axis')
% hold on
% yyaxis("left")
% ylabel("MAPE")
% plot(timeline/3600,MAPE(1,:))


% Seasonal analysis 
% It is evident from figures 2 and 3 that the error function is time
% dependent. However, it is not clear that this is the only factor at play.
% Namely, because the shape of the swarm changes as a funciton of time
% also, it may be that the error function is swarm-density dependent
% instead, or more strongly so. The following analyisis is carried out to
% verify if a correlation exists between swarm density and the Mean Average
% Percentage Error

% 1. Get the swarm density as a function of time  
% get ADJ in cylindrical coordinates (in-line density is the most
% representatve of swarm structure)
ADJcyl = cell(1,length(timeline));
for t=1:length(timeline)
    ADJcyl{t} = getADJ(dataCyl,t);
end
% Establish swarm boundaries at each instant. 
swarmScale = zeros(3,length(timeline)); 
for t=1:length(timeline)
    swarmScale(1,t) = max(max(ADJcyl{t}{1}));%Theta
    swarmScale(2,t) = max(max(ADJcyl{t}{2}));%Rho
    swarmScale(3,t) = max(max(ADJcyl{t}{3}));%z
end
% Get the average density as a function of time (n/envelope)
n = height(Data{1});
swarmDensity = n./swarmScale; 
%% Plot swarm density with MAPE
figure(3)
hold on
title('Correlation betwen MAPE and density along z axis')
hold on
yyaxis("left")
ylabel("MAPE")
plot(timeline,MAPE(1,:))
yyaxis("right")
plot(timeline,swarmDensity(3,:))
%legend('MAPE','density')
xlabel("time [h]")
ylabel("density [n/m]")
%xlim([0,timeline(length(timeline))])
%% draw period lines
period = 6328/3600;
for Tau=1:5 
    xline(Tau*period)
    xline(Tau*0.5*period)
end
% Get the MAPE correlation for each of the other two axis 
figure(4)
hold on
title("Correlation of MAPE to density along \Delta \Theta")
yyaxis("left")
ylabel("MAPE")
plot(timeline,MAPE(1,:))
yyaxis("right")
ylabel("density n/rad")
plot(timeline,swarmDensity(1,:));
legend('MAPE','density')
xlabel("time [h]")
xlim([0,timeline(length(timeline))])

figure(5)
hold on
title("Correlation of MAPE to density along \Delta \rho")
yyaxis("left")
ylabel("MAPE")
plot(timeline,MAPE(1,:))
yyaxis("right")
ylabel("density [n/m]")
plot(timeline,swarmDensity(2,:));
legend('MAPE','density')
xlabel("time [h]")
xlim([0,timeline(length(timeline))])

clear n t ADJcyl adjCyl


%% -------------------------FUNCTIONS-------------------------------- %% 
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
    %disp("WARNING: x2 is not an integer, so p has been recalculated")
    %disp(p)       
end
x1 = n-x2;
%fprintf('X1 has %i elements \n',x1)
%fprintf('X2 has %i elements \n',x2)

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
    % generate a group size "c" of at least two satellites, and at most
    % half the list size. 
    if cCount == 0
        %Block group size 1 and anything greater than half the size of the
        %test set
        allowedC = setdiff(1:length(allowedList),[1 x2/2:x2]);
        %Get group size
        c = allowedC(randi(length(allowedC))); 
    else
        allowedC = setdiff(1:length(allowedList),[1,length(allowedList)-1]);
        %don't allow c to be greater than number of slots left
        c = allowedC(randi(length(allowedC))); 
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

% GPPERFORMANCE: is a functionthat takes a data set D, a test set pointer
% vector, the set of y predictions at the locations pointed to by the test
% set, and the time during which the observerations where made as
% arguments. The function returns the MAE, MAPE, and RMSE for those
% predictions given the set of known measurements associated with the test
% set. 
function [MAE,MAPE,RMSE] = gpPerformance(D,testSet,y_pred_1,T)
    % Calculate stats
    % Root Mean Square Error (RMSE) for the test set
    RMSE = sqrt((1/length(testSet))*sum((y_pred_1 - D{T}(testSet,7)).^2));
    % Mean Absolute Error
    MAE = (1/length(testSet))*sum(abs(y_pred_1 - D{T}(testSet,7)));
    % Mean Aboslute Percentage Error
    MAPE = (1/length(testSet))*sum(abs((y_pred_1 - D{T}(testSet,7))./D{T}(testSet,7)));
end

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


