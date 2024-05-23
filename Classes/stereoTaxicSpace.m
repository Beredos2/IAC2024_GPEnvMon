%---------------------------StereoTaxicSpace------------------------------%
%A stereotaxic space is a normalized space created to index locations in
%real space. 
%Properties
%dl = side length of a cartesian cube (cell side lenght) 
%Methods
%getIndex(R,xbias,ybias,zbias,spaceObj): returns the cell index. Bias is used to
%floor the function up or down along each axis. 
classdef stereoTaxicSpace < handle
    properties
        dl
    end
    methods
        %CONSTRUCTOR
        function spaceObj = stereoTaxicSpace(dl)
            spaceObj.dl = dl; %Input dependent units
        end
        %GET INDEX: is a function that returns the index of a position R
        %within the stereotaxic space
        function index = getIndex(R,xbias,ybias,zbias,STS)
            %Modulate
            index = R./STS.dl; 
            %Identify singularities
            if index(1)==0
                if xbias >= 0
                    index(1) = 1;
                else
                    index(1) = -1;
                end
            end
            if index(2) == 0
                if ybias >= 0
                    index(2) = 1; 
                else
                    index(2) = -1;
                end
            end 
            if index(3) == 0
                if zbias >= 0
                    index(3) = 1;
                else
                    index(3) = -1;
                end
            end         
            
            %Check sign and eliminate decimals
            if index(1) > 0
                index(1) = ceil(index(1)); 
            else
                index(1) = floor(index(1)); 
            end 
            
            if index(2) > 0 
                index(2) = ceil(index(2));
            else 
                index(2) = floor(index(2)); 
            end

            if index(3) > 0
                index(3) = ceil(index(3));
            else 
                index(3) = floor(index(3));
            end
        end
        %CHECK r
        function trueFalse = checkR(R,index,STS)
            %create boundaries
            if index(1)>0
                xmax = index(1)*STS.dl;
                xmin = index(1)*STS.dl-STS.dl; 
            else
                xmax = index(1)*STS.dl+STS.dl;
                xmin = index(1)*STS.dl;
            end
            if index(2)>0
                ymax = index(2)*STS.dl;
                ymin = index(2)*STS.dl-STS.dl;
            else
                ymax = index(2)*STS.dl+STS.dl;
                ymin = index(2)*STS.dl;
            end
            if index(3)>0
                zmax = index(3)*STS.dl;
                zmin = index(3)*STS.dl-STS.dl; 
            else
                zmax = index(3)*STS.dl + STS.dl;
                zmin = index(3)*STS.dl; 
            end
            
            % X = [xmin xmax]
            % Y = [ymin ymax]
            % Z = [zmin zmax]

            %check R
            if R(1) <= xmax && R(1) >= xmin
                trueFalse_X = 1;
            else 
                trueFalse_X = 0;
            end
            if R(2) <= ymax && R(2) >= ymin
                trueFalse_Y = 1;
            else 
                trueFalse_Y = 0;
            end
            if R(3) <= zmax && R(3) >= zmin
                trueFalse_Z = 1;
            else 
                trueFalse_Z = 0;
            end

            if trueFalse_X == 0 || trueFalse_Y == 0 || trueFalse_Z == 0
                trueFalse = 0; 
            else
                trueFalse = 1; 
            end
            %[trueFalse_X,trueFalse_Y,trueFalse_Z]
        end
        %WhatCell: is a function that takes a position and an STS object,
        %and returns the cell index for that position (always biased up)
        % function index = whatCell(R,STS)
        %     index = [0 0 0];
        %     %R is in cell (x,y,z), when R_z<x, R_y<y, and R_z<z
        %     %Normalize R
        %     R = R/STS.dl; 
        %     %Check for singularities 
        %     if R(1)==0
        %         R(1) = 1; 
        %         disp("singularity in x")
        %     end
        %     if R(2)==0
        %         R(2)= 1;
        %         disp("singularity in y")
        %     end
        %     if R(3) == 0
        %         R(3) = 1; 
        %         disp("singularity in z")
        %     end
        %     %Find upper limmits
        %     index(1) = ceil(R(1));
        %     index(2) = ceil(R(2)); 
        %     index(3) = ceil(R(3)); 
        % end

    end
end