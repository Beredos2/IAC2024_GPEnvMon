function ADJ = getADJ(myData,time)
%GETADJ is a function that takes a data set of the form
%myData{timeIndex}(id,[x,y,z][M_x,M_y,M_z][M]) and time, and returns the weighted
%adjacency matrices for the spatial component of the set for that time, were the weight of
%each adjacency variable is the distance bewteen two points of measurement.
if nargin == 2
    %Create a buffer to be filled out for each ADJ
    ADJ_theta = zeros(height(myData{time})); 
    ADJ_rho = zeros(height(myData{time})); 
    ADJ_z = zeros(height(myData{time})); 
    
    %Fill in weights
    for i=1:height(myData{time})    
        for j = i:height(myData{time})
            %THETA
            ADJ_theta(i,j) = myData{time}(j,1)-myData{time}(i,1);
            ADJ_theta(j,i) = -ADJ_theta(i,j); 
            %RHO
            ADJ_rho(i,j) = myData{time}(j,2) - myData{time}(i,2);
            ADJ_rho(j,i) = -ADJ_rho(i,j);
            %Z
            ADJ_z(i,j) = myData{time}(j,3) - myData{time}(i,3); 
            ADJ_z(j,i) = -ADJ_z(i,j);
        end
    end
    
    ADJ = cell(1,3);
    ADJ{1} = ADJ_theta;
    ADJ{2} = ADJ_rho; 
    ADJ{3} = ADJ_z; 
    
    %disp("your adjacency matrices are arranged in a cell structure such that {ADJ_theta,ADJ_rho,ADJ_z}")
elseif nargin == 1
    %Create a buffer to be filled out for each ADJ
    ADJ_theta = zeros(height(myData)); 
    ADJ_rho = zeros(height(myData)); 
    ADJ_z = zeros(height(myData)); 
    
    %Fill in weights
    for i=1:height(myData)    
        for j = i:height(myData)
            %THETA
            ADJ_theta(i,j) = myData(j,1)-myData(i,1);
            ADJ_theta(j,i) = -ADJ_theta(i,j); 
            %RHO
            ADJ_rho(i,j) = myData(j,2) - myData(i,2);
            ADJ_rho(j,i) = -ADJ_rho(i,j);
            %Z
            ADJ_z(i,j) = myData(j,3) - myData(i,3); 
            ADJ_z(j,i) = -ADJ_z(i,j);
        end
    end
    
    ADJ = cell(1,3);
    ADJ{1} = ADJ_theta;
    ADJ{2} = ADJ_rho; 
    ADJ{3} = ADJ_z; 
    
    %disp("your adjacency matrices are arranged in a cell structure such that {ADJ_theta,ADJ_rho,ADJ_z}")
end
