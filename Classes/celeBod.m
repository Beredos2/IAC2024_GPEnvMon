%----------------------CELESTIAL BODY-------------------------------------%
%celeBod(R,V,Theta,w) instantiates a celestial body in location R with
%respect to the frame of reference, velocity V, orientation Theta, rotation
%rate w, where R,V,Theta, and w are vectors in three dimensions, and whose
%values are measured with respect to an inertial frame of reference.   

classdef celeBod < handle
    properties
        R       %Position m
        V       %Velocity m/s
        Theta_0   %Orientation rad
        Theta_t %Orientation history (array)
        w       %Angular rate rad/s
        radius  %Planetary radius in km
    end

    methods
        %------------------------CONSTRUCTOR-------------------------------
        function CelestialBody = celeBod(R,V,Theta,w,radius,C_obj)
            CelestialBody.R = R;
            CelestialBody.V = V;
            CelestialBody.Theta_0 = Theta;
            CelestialBody.w = w;
            CelestialBody.radius = radius;
            %Calculate orientation history
            CelestialBody.Theta_t = zeros(length(C_obj.timeline),3); %Orientations stored in a nTx3 array, where the columns are Theta_X,Theta_Y,Theta_Z in the ECI frame, and states throughout time are read top to bottom.
            for t=1:length(C_obj.timeline)
                CelestialBody.Theta_t(t,:) = getOrientations(CelestialBody,C_obj.timeline(t));
            end
                
        end
        %-------------------------DYNAMICS---------------------------------
        function Theta_t = getOrientations(CelestialBody,t) %given in radians
            Theta_x = CelestialBody.Theta_0(1)+CelestialBody.w(1)*t;
            Theta_y = CelestialBody.Theta_0(2)+CelestialBody.w(2)*t;
            Theta_z = CelestialBody.Theta_0(3)+CelestialBody.w(3)*t;
            Theta_t = [Theta_x,Theta_y,Theta_z];
        end
        %------------------------MANETIC FIELD----------------------------%
        function [XYZ,F] = planetMagField(P_obj,height,latitude,longitude,epoch)
            % epoch to decimal year
            d = datetime(epoch,'ConvertFrom','epochtime','Format','yy/MM/dd');
            d = decyear(d);
            p = 1; 
            % Check for an event
            if p == 0%rand>p
                %call m=wrldmagm
                [XYZ,~,~,~,F] = wrldmagm(height,latitude,longitude,d);
            else
                %call wrldmagm*
                [XYZ,F] = myWrldmagm(P_obj,height,latitude,longitude,d,epoch);
            end             
        end

        %------------------------myWrldmagm----------------------------%
        %returns the world magnetic model prediction of the magnetic field
        %vector and intensity, with an added signal. 
        function [XYZ,F] = myWrldmagm(P_obj,height,latitude,longitude,d,epoch)
            %Get the baseline magnetic field
            [XYZ] = wrldmagm(height,latitude,longitude,d);
            %Define signal constants 
            A1 = 0.333/2; %x axis
            A2 = 0.455/2; %y axis
            A3 = 3.66/2;  %z axis
            A = [A1,A2,A3]; 
            clear A1 A2 A3
            freq1 = 1.25*pi/1800; %Every hour and a quarter. 
            freq2 = 5*pi/30; %5 minutes
            freq3 = 10*2*pi; %10 seconds
            freqL = 10e5; 
            freq = [freq1,freq2,freq3]; 
            clear freq1 freq2 freq3
            %Calculate the contribution to each component
            Xstar = XYZ(1)*A(1)*(cos(freq(1)*epoch)+cos(freqL*latitude));
            Ystar = XYZ(2)*A(2)*(cos(freq(2)*epoch)+cos(freqL*latitude));
            Zstar = XYZ(3)*A(3)*(cos(freq(3)*epoch)+cos(freqL*latitude));
            %Add Mstar to M
            XYZ(1) = XYZ(1) + Xstar; 
            XYZ(2) = XYZ(2) + Ystar;
            XYZ(3) = XYZ(3) + Zstar;
            F = norm(XYZ); 
        end        
    end
end 