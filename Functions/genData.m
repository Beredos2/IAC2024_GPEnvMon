%genData(drone,clock) is a function that takes a drone and clock objects as
%arguments, and returns an array that contains the simulated data
%collection for that drone throughout its lifecycle [P_ECEF,|M|,t_dect],
%where, 
%P_ECEF are the planetary coordinates (longitude,latitude,R) of the
%drone
%|M| is the magnetic field intensity 
%t_dect is the time of detection 

function genData(D_obj,C_obj)
    %Get space coordinates
    R_t = D_obj.r; 
    %Get time coordinates
    T = C_obj.timeLine;
    %Convert ECI space coordinates to ECEF coordinates
    %1.Find the Earth's roation when each measurement was taken
    %2.Construct a rotation matrix
    %3.Rotate all timeline vectors 


end