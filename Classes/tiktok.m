%Clock
%This class is used to store properties related to the simulation duration,
%step size, and duration calculation. 

classdef tiktok < handle
    properties 
        epoch_0     %Mission start  (epoch in seconds)
        epoch_f     %End of mission (epoch in seconds)
        dT          %stepsize in s
        nT          %total number of steps
        T           %mission duration
        timeline    %vector that stores discreet time measurements throughout mission lifecycle
    end

    methods
        %CONSTRUCTOR
        function C_obj = tiktok(epoch_0,epoch_f,dT)
            C_obj.epoch_0 = epoch_0;
            C_obj.epoch_f = epoch_f;
            C_obj.T = epoch_f-epoch_0;           
            C_obj.nT = ceil(C_obj.T/dT); %calculate number of intervals needed to the nearest higher integer
            C_obj.T = C_obj.nT*dT; %Recalculate mission duration
            C_obj.dT = dT;
            C_obj.epoch_f = epoch_0+C_obj.T; %assign new final epoch based on updated duration
            fprintf('To keep whole units of time, mission duration has been update to %d, so that the mission end epoch was updated to %d \n',C_obj.T,C_obj.epoch_f)
            %Construct timeline
            C_obj.timeline = zeros(C_obj.nT,2); %The timeline contains both the 0 to T time in timeline(:,1), used for mission reference, and the epoch values throughout the mission duration in timeline(:,2)
            C_obj.timeline(:,1) = linspace(0,C_obj.T,C_obj.nT);
            %Construct epoch timeline
            C_obj.timeline(:,2) = linspace(C_obj.epoch_0,C_obj.epoch_f,C_obj.nT);
        end
    end
   
end
