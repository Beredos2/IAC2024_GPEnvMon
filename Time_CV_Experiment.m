%Run Data revovery and 
%GP_1 before running this script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBEXTRACTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- You can adapt the sub exctraction

Y=zeros(8150,1);
Z=zeros(8150,1);
X1=zeros(8150,1);
X2=zeros(8150,1);
X3=zeros(8150,1);

% Change timeline to hours
for i = 1:50
start=(i*163)-162;
ends=i*163;
%extract time wise measurements
% extract spatial coordinates
y= [Data{i}(:,7)];   %magnetic field intensity
x1=[Data{i}(:,1)];
x2=[Data{i}(:,2)];                   
x3= [Data{i}(:,3)];  % 3d-coordinate-x

%create vectors of spatial coordinates
X1(start:ends)=x1;
X2(start:ends)=x2;
X3(start:ends)=x3;
Y(start:ends)=y;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT OF THE SUBEXTRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot magnetic field intensity as a function of Z 
figure(1)
plot(X1,Y,'bo')
%plot(1:163,Y,'bo')
%legend('predicted values','measured values')
xlabel('X1}')
ylabel('Magnetic field intensity (nT)')
%
figure(2)
plot(X2,Y,'bo')
%plot(1:163,Y,'bo')
%legend('predicted values','measured values')
xlabel('X2')
ylabel('Magnetic field intensity (nT)')


figure(2);
scatter3(X1, X2, X3, 'filled');
title('3D Scatter Plot');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROSS VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_T=(2:10:50);
%time point t=2,t=10,t=30,t=40,t=50
num_of_time_wise_exp=10;
ExperimentCV_1D=zeros(length(sample_T),num_of_time_wise_exp);
ExperimentCV_2D=zeros(length(sample_T),num_of_time_wise_exp);
ExperimentCV_3D=zeros(length(sample_T),num_of_time_wise_exp);


for t = 1:5
    
t1=sample_T(t);
t_index_start=(t1*163)-162;
t_index_end=t1*163;
MAE_1D=zeros(num_of_time_wise_exp,1);
MAE_3D=zeros(num_of_time_wise_exp,1);
MAE_3D=zeros(num_of_time_wise_exp,1);
for rep = 1:num_of_time_wise_exp
% 
y_t=Y(t_index_start:t_index_end);
x1_t=X1(t_index_start:t_index_end);
x2_t=X2(t_index_start:t_index_end);
x3_t=X3(t_index_start:t_index_end);
%Create test and training set
%data_set=size, missing_fraction,averge_missing_lentgth, Boundary
test_idx=BB_idx(163,0.2,5,10);
train_idx= setdiff(1:163', test_idx)';
%predictors matrix
X=[x1_t,x2_t,x3_t];


Model_1D_Z = fitrgp(X(train_idx,3), y_t(train_idx),'KernelFunction','squaredexponential');
Model_2D_X_Y = fitrgp(X(train_idx,1:2), y_t(train_idx),'KernelFunction','squaredexponential');
Model_3D = fitrgp(X(train_idx,1:3), y_t(train_idx),'KernelFunction','squaredexponential');

[yPred_1D, ySD_1D, yInt_1D] = predict(Model_1D_Z, X(test_idx,3));
[yPred_2D_X_Y, ySD_2D_X_Y, yInt_2D_X_Y] = predict(Model_2D_X_Y, X(test_idx,1:2));
[yPred_3D , ySD_3D , yInt_3D ] = predict(Model_3D , X(test_idx,1:3));

MAE_1D(rep)=mean(abs(yPred_1D-y_t(test_idx)))
MAE_2D(rep)=mean(abs(yPred_2D_X_Y-y_t(test_idx)))
MAE_3D(rep)=mean(abs(yPred_3D-    y_t(test_idx)))

end

ExperimentCV_1D(t,1:10)=MAE_1D;
ExperimentCV_2D(t,1:10)=MAE_2D;
ExperimentCV_3D(t,1:10)=MAE_3D;
end

figure(1)
plot(X(train_idx,3), y_t(train_idx),'bo')
hold on 
scatter(X(test_idx,3),y_t(test_idx),'*')
plot(X(test_idx,3),yPred_1D,'ro')
legend('train values','test values','prediction')
xlabel('X3')
ylabel('Magnetic field intensity (nT)')



