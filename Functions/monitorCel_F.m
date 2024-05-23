%[T F] = monitorCel(C_obj,S_network,STS_obj,index): is a function that takes a clock
%object, a swarm's network, a steretaxic space object, and an index, and
%returns the average magnetic field intensisty of the indexed cell measured
%for all time. 

function [T, F] = monitorCel_F(C_obj,S_network,STS_obj,index)
%generate a storage matrix with the dimensions of time and number of drones
T = C_obj.timeline(:,2);
M = ones(length(T),length(S_network)); 
M = -M; 
%Fill out storage matrix only when measurements are taken in the indexed
%cell
for i = 1:length(S_network)
    for j = 1:length(T)
        bool = checkR(S_network{i}.r(:,j),index,STS_obj);
        if bool == 1
            M(j,i) = S_network{i}.D_m(j,4);
        end
    end
end
%disp(M(1:10,:))
%Get cell average for each time step 
F = zeros(height(M),1);
for i = 1:height(M)
    x = 0; 
    for j = 1:width(M)
        if M(i,j)>0
            F(i) = F(i)+M(i,j);
            %disp(M(i,j))
            x = x+1;
        end
        if x > 0
            F(i) = F(i)/x; 
        else
            fprintf("no measurements for this cell at %d \n",T(i))
        end
    end
end

