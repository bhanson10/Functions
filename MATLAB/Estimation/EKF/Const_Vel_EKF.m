clear all; close all; clc;

% Const_Vel_EKF.m
% Ben Hanson
% 1/5/2024 
%
% This example is pulled directly from the MATLAB website, but uses the
% formulae of the EKF rather than the built-in "trackingEKF" function


rng(2022); % For repeatable results
dt = 0.5;  % Time step 
T  = 20;    % Period of propagation
x0 = [30; 40; 1; 1];        % [x;y;vx;vy]
Q  = diag([0; 0; .01; .01]); % Process noise matrix
R  = diag([2e-6;1]);         % Measurement noise matrix. Units are m^2 and rad^2.

% True States and Measurements
[x, z, tx, tz] = get_truth_measurements('x0', x0, ...
                                        'dt', dt, ...
                                        't', [dt,T], ...
                                        'f', @f, ...
                                        'h', @h, ...
                                        'Q', Q, ...
                                        'R', R, ...
                                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation

% Extended Kalman Filter
xest0 = [32; 38; 0; 0]; Pest0 = diag([2,2,1e3,1e3]);
[xest,Pest,tspan] = EKF('xest0', xest0, ...
                        'Pest0', Pest0, ...
                        'dt', dt, ...
                        'z', {tz,z}, ...
                        'Q', Q, ...
                        'R', R, ...
                        'f', @f, ...
                        'h', @h, ...
                        'F', @F, ...
                        'H', @H, ...
                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation/EKF

figure(1)
subplot(2,1,1)
plot(tz,z(1,:)*180/pi,'LineWidth',1)
xlabel("time (s)")
ylabel("angle (deg)")
title("Angle and Range")
subplot(2,1,2)
plot(tz,z(2,:),'LineWidth',1)
xlabel("time (s)")
ylabel("range (m)")

figure(2); hold on; axis equal;  
cad = round(linspace(1,length(tx),6));
for i=cad
    Pest1 = Pest{i};
    plot_gaussian_ellipsoid(xest(1:2,i),Pest1(1:2,1:2),3);
end
plot(x(1,:),x(2,:),"r-",'LineWidth',1,"DisplayName","True Trajectory");
plot(xest(1,:),xest(2,:),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
scatter(x(1,cad),x(2,cad),'r', 'filled','HandleVisibility','off');
scatter(xest(1,cad),xest(2,cad),'b', 'filled','HandleVisibility','off');
xlim([25,60])
ylim([30,65])
xlabel("x (m)")
ylabel("y (m)")
legend("Location","northwest")
title("True Trajectory vs Estimated Trajectory")
axis square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(x,const)
    x1 = [x(3); x(4); 0; 0]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x)
    z = [atan(x(2)/x(1)); norm([x(1) x(2)])]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = F(x,const)
    Fx = [0 0 1 0;         
          0 0 0 1;
          0 0 0 0;
          0 0 0 0];
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = H(x)
    Hx = [-x(2)/(x(1)^2+x(2)^2)    x(1)/(x(1)^2+x(2)^2)     0 0;
          x(1)/sqrt(x(1)^2+x(2)^2) x(2)/sqrt(x(1)^2+x(2)^2) 0 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%