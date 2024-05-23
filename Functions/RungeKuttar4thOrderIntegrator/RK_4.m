%==========================================================================
%  RK4 WITH POSITION AND VELOCITY
%==========================================================================
%> This function numerically solve equation of orbital motion using 
%> the fourth order Runge-Kutta method.
%>
%> @author Ozan Kilic
%> Middle East Technical University, Space Geodesy Division
%> contact: kilic.ozan@metu.edu.tr
%==========================================================================
function [X_RK] = RK_4(X,h,steps)

X_init = X;

for t=1:steps
    
X=rk4orbit(t,X,h);
XX(t)=X(1);
YY(t)=X(2);
ZZ(t)=X(3);
VX(t)=X(4);
VY(t)=X(5);
VZ(t)=X(6);

 % plot3(XX,YY,ZZ);   % Plot numerical orbit
 % pause(0.0000000005)

end

% title('Simulated Orbit in ECI')

% Create PV Matrix
X_RK = vertcat(XX,YY,ZZ,VX,VY,VZ);
X_RK = horzcat(X_init,X_RK);

end


function X = rk4orbit(t,X,h)
% Numerical Solution, Runge-Kutta 4th Order
    k1 = forbit(t,X);
    k2 = forbit(t+h/2,X+k1*h/2);
    k3 = forbit(t+h/2,X+k2*h/2);
    k4 = forbit(t+h,X+k3*h);

% Step forward in time
X = X+(h/6)*(k1+2*k2+2*k3+k4);

end

function Xdot = forbit(t,X)

    mu = 3.986004418*10^14;    % Gravitational parameter of Earth (m^3/s^2)
    r = X(1:3);                % Position (m)
    v = X(4:6);                % Velocity (ms^2)

    dr = v;
    dv = (-mu/(norm(r))^3).*r; % Newton's law of gravity
    Xdot = [dr; dv];
    
end