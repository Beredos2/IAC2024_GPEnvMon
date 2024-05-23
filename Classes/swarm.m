%--------------------------SWARM------------------------------------------%
%A swarm is an object that is composed of two or more drones. It inherits
%properties from the handle class. 
%
%CONSTRUCTOR
%The class constructor takes 4 inputs (ID,type,systems,IC)
%1.ID: is a number that uniquely identifies the swarm
%2.type: is a string that describes the function of the swarm
%3.systems: is a table that describes the size and requirements of each
%   sub-swarm. The table headers are [dTypes,dPop,dNeed], where they are the
%   type of drone, the drone population, and the required drones for mission
%   success respectively. 
%4.IC: are the initial conditions of the swarm. These are part of a vector
%   that stores the starting keplerian coordinates [a,e,i,Omega,w,M], where
%   -a: semi-major axis of the deployer's orbit [km]
%   -e: eccentricity of the deployer's orbit
%   -i: inclination of the deployer's orbit [rad]
%   -Omega: right ascension of the ascending node """ [rad]
%   -w: angle of periapse [rad]

classdef swarm < handle
    properties
    %COMPONENTS 
    r_reff     %is a virtual satellite that follows an unperturbed orbit from deployment and serves as a refference for all other members of the swarm
    network    %is an array of objects containing all "real" drones
    %DEFINITION
    ID         %identification digit
    type       %swarm type
    systems    %stores a table that describes each system based on dTypes, dPop, and dNeed
    Kepler     %stores the initial conditions of the swarm at deployment: [a,e,i,Omega,w,M]
    t_0        %stores the epoch at time of deployment
    %STATE
    %adjMatrix  %stores the range relationships between satellites :
    %consider exporting adjacency matrix to data anlysis rather than in
    %object
    age        %measures time since deployment
    %epsilon_g  %measures the gathering efficiency of the swarm
    end
    
    methods
        %CONSTRUCTOR: receives ID,type, systems, and ICs as inpuTs, and
        %generates a swarm with those properties. 
        function S_obj = swarm(ID,type,systems,Kepler,t_0)
            %Assign inputs to swarm properties
            S_obj.ID = ID; 
            S_obj.type = type;
            S_obj.systems = systems;
            S_obj.Kepler = Kepler;
            S_obj.t_0 = t_0;
            
            buildNetwork(S_obj); %Assemble swarm
        end

        %BUILD NETWORK: instantiate all the members of a sawrm, and assigns
        %them a place in the network. The drone arrays are arranged 
        %horizontally on the basis of type, and then stacked vertically. 
        %Because the sub-swarm sizes may vary, the vertical arrays are
        %nested in the horizotnally arranged cells, so that two cell calls
        %are required to reach any drone. 
        function buildNetwork(S_obj) 
            %1. Create a cell array capable of storing the network
            S_network = cell(1,height(S_obj.systems)); 
            %2. Loop through each node type
            for i=1:height(S_obj.systems)%j indicates the drone types
                %2.1 Create a cell array to store all drones of each type
                snetwork = cell(S_obj.systems.dPop(i),1);
                %2.2 Assign drones to the network with unique identifiers
                for j=1:S_obj.systems.dPop(i) %where j indicates the drone in a sub-swarm
                    snetwork{j} = drone([i, j],S_obj.systems.dType(i),S_obj.Kepler);
                end
                %2.3 Assign sub network to the general swarm
                S_network{i} = snetwork;
            end
            S_obj.network = S_network;     
        end
        
        %DEPLOY: deployment applies a delta vee to each drone in the
        %netowrk, causing its initial conditions to change before
        %propagation. NOTE: take into account that multiple stages need to
        %be propagated separately, since each propagation occurs throughout
        %the entire time of the clock. 
        function deployBrood(S_obj) %NOTEL: as written, the function generates the velocities internally. While this is easier than providing them as an imput, it doesn't allow the user to tune the velocity distribution. Future update will change this
                %1. Generate a set of random velocities
                rng(0,'twister'); %initializes random number generator so that results are reproducible
                dv = cell(height(S_obj.systems)); %buffer to store velocities
                
                for i=1:length(dv)
                    dV_x = normrnd(0.1,0.001,[1,S_obj.systems.dPop(i)]); %km/s
                    dV_y = normrnd(0.1,0.001,[1,S_obj.systems.dPop(i)]); %km/s
                    dV_z = normrnd(0,0.0001,[1,S_obj.systems.dPop(i)]); %km/s
                    dv{i} = [dV_x;dV_y;dV_z]; 
                end
                %2. Add velocity to the V_0 of each drone
                for i=1:length(dv)
                    for j=1:S_obj.systems.dPop(i)
                        S_obj.network{i}{j}.v(1) = S_obj.network{i}{j}.v(1)+dv{i}(1,j);
                        S_obj.network{i}{j}.v(2) = S_obj.network{i}{j}.v(2)+dv{i}(2,j);
                        S_obj.network{i}{j}.v(3) = S_obj.network{i}{j}.v(3)+dv{i}(3,j);
                    end
                end
        end
        
        %getAdjMatrix_R: is a method that returns the adjacency matrix
        %describing the position of drones at a given time, for a given
        %swarm. 
        function adjM_R = getAdjMatrix_R(S_obj,C_obj,T)
            %Get model based time
            i=1;
            while T>C_obj.timeline(i,1)
            i=i+1;
                if i>height(C_obj.timeline)
                    disp("time is not contained in timeline")
                    break
                end
            end
            T = C_obj.timeline(i-1,1);
            fprintf("time changed to %d to fit model time records",T)
            %Fill the adjacency matrix for each network
            A = cell(1,length(S_obj.network));
            for i=1:length(S_obj.network)
                A{i} = zeros(length(S_obj.network{i}));
                for j=1:height(A) %reference counter
                    for k=(1+j):height(A) %referred counter
                        A{i}(j,k) = norm(S_obj.network{i}{k}.r(:,T)-S_obj.network{i}{j}.r(:,T));
                        A{i}(k,j) = A{i}(j,k); %symmetry
                    end
                end
            end
            adjM_R = A; 
        end
    end

end

    
