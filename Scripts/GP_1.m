% GP_01: In this script we use the fitGP function to predict the value of
% the mangetic field intensity given any set of spatial coordiantes. For
% the training data, we recover data from an MDIPS demonstration
% simulation by running DataRecovery_PolarCoordinates script.
clear 
clc
addpath("C:\Users\2700783t\OneDrive - University of Glasgow\InMotu\MATLAB\Gaussian Process Regression")
DataRecovery

% Change timeline to hours
timeline = timeline/3600; 
tindex = 2000; 

% fit a GP model for timeline index 2000  (correlate |M| to Z
X = [Data{tindex}(:,3)]; %position
Y = Data{tindex}(:,7);     %magnetic field intensity

Model_1 = fitrgp(X,Y,'KernelFunction','squaredexponential');

y_pred_1 = resubPredict(Model_1); 

% fit a GP model for timeline index 2000  (correlate |M| to X
X = [Data{tindex}(:,1)]; %position
Y = Data{tindex}(:,7);     %magnetic field intensity

Model_2 = fitrgp(X,Y,'KernelFunction','squaredexponential');

y_pred_2 = resubPredict(Model_2); 



%% 


%% Plot the data in R^3 space 
figure(10) 
scatter3(Data{tindex}(:,1),Data{tindex}(:,2),Data{tindex}(:,3),10,Data{tindex}(:,7),'filled')
colorbar
set(gca,'colorscale','log')
xlabel('X_{ECI} [m]')
ylabel('Y_{ECI} [m]')
zlabel('Z_{ECI} [m]')
title('Magnetic field magnitude as a function of position')

%% Plot magnetic field intensity as a function of Z 
figure(2)
plot(Data{tindex}(:,3),y_pred_1,'r.')
hold on
plot(Data{tindex}(:,3),Y,'bo')
%plot(1:163,Y,'bo')
legend('predicted values','measured values')
xlabel('Z_{ECI}')
ylabel('Magnetic field intensity (nT)')



