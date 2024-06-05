clear all; close all; clc;

% Lorenz_EKF.m
% Ben Hanson, 2024
%
% Applying the EKF to the Lorenz system

% Initial conditions
initialize_figures(); 
rng(2024); % Random Number Generator
const.sigma=4; const.b=1; const.r=48; % Lorenz Parameters
T   = 1; % Period of propagation
dt  = 0.02; % Time step or [True time step, Measurement time step]
x0  = [-11.5; -10; 9.5]; % [x;y;z]
Q   = diag([.01; .01; .01;]); % Process noise matrix
R   = diag([1; 0.1; 0.1]); % Measurement noise matrix

% True States and Measurements
[x, z, tx, tz] = get_truth_measurements('x0', x0, ...
                                        'dt', dt, ...
                                        't', [dt,T], ...
                                        'f', @f, ...
                                        'h', @h, ...
                                        'Q', Q, ...
                                        'R', R, ...
                                        'const', const, ...
                                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation

% Extended Kalman Filter - Explicit Euler
[xest,Pest,txest] = EKF('xest0', x0, ...
                        'Pest0', diag([1,1,1]), ...
                        'dt', dt, ...
                        'z', {tz, z}, ...
                        'Q', Q, ...
                        'R', R, ...
                        'f', @f, ...
                        'h', @h, ...
                        'F', @F, ...
                        'H', @H, ...
                        'const', const, ...
                        'method', 'EE', ...
                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation/EKF

% Extended Kalman Filter
[xest2,Pest2,txest] = EKF('xest0', x0, ...
                          'Pest0', diag([1,1,1]), ...
                          'dt', dt, ...
                          'z', {tz, z}, ...
                          'Q', Q, ...
                          'R', R, ...
                          'f', @f, ...
                          'h', @h, ...
                          'H', @H, ...
                          'const', const, ...
                          'method', 'RK4', ...
                          'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation/EKF


% Plotting Background behavior
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) f(Y,const), [0,50], x0, options);
figure(2); 
plot3(Y(:,1),Y(:,2),Y(:,3),'Color', '[0 1 0 0.5]','linewidth',0.3, 'HandleVisibility','off');
drawnow;

figure(1)
subplot(3,1,1)
plot(tz,z(1,:),'LineWidth',1)
xlabel("time (s)")
ylabel("range (m)")
subplot(3,1,2)
plot(tz,z(2,:)*180/pi,'LineWidth',1)
xlabel("time (s)")
ylabel("\theta (deg)")
subplot(3,1,3)
plot(tz,z(3,:)*180/pi,'LineWidth',1)
xlabel("time (s)")
ylabel("\phi (deg)")

figure(2); hold on;
cad = round(linspace(1,length(tx),6));
plot3(x(1,:),x(2,:),x(3,:),"k-",'LineWidth',1,"DisplayName","True Trajectory");
scatter3(x(1,cad),x(2,cad),x(3,cad),'k', 'filled','HandleVisibility','off');

plot3(xest(1,:),xest(2,:),xest(3,:),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory, EE");
scatter3(xest(1,cad),xest(2,cad),xest(3,cad),'b', 'filled','HandleVisibility','off');
plot3(xest2(1,:),xest2(2,:),xest2(3,:),"r--",'LineWidth',1,"DisplayName","Estimated Trajectory, RK4");
scatter3(xest2(1,cad),xest2(2,cad),xest2(3,cad),'r', 'filled','HandleVisibility','off');

for i=round(linspace(1,length(tx),6))
    p.color = [0 0 1];
    plot_gaussian_ellipsoid(xest(1:3,i),Pest{i},3,p);
    p.color = [1 0 0];
    plot_gaussian_ellipsoid(xest2(1:3,i),Pest2{i},3,p);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=f(x,const)                          
    f=[const.sigma*(x(2)-x(1));  -x(2)-x(1)*x(3);  -const.b*x(3)+x(1)*x(2)-const.b*const.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x)
    z = [sqrt(x(1)^2+x(2)^2+x(3)^2); atan2(x(2),x(1)); acos(x(3)/sqrt(x(1)^2+x(2)^2+x(3)^2))]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = F(x,const)
    Fx = [-const.sigma const.sigma 0       ;         
          -x(3)        -1          -x(1)   ;
          x(2)         x(1)        -const.b];
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = H(x)
    Hx = [x(1)/sqrt(x(1)^2+x(2)^2+x(3)^2) x(2)/sqrt(x(1)^2+x(2)^2+x(3)^2) x(3)/sqrt(x(1)^2+x(2)^2+x(3)^2);
          -x(2)/(x(1)^2+x(2)^2)           x(1)/(x(1)^2+x(2)^2)            0;
          x(1)*x(3)/((x(1)^2+x(2)^2+x(3)^2)^(3/2)*sqrt(1-(x(3)^2/(x(1)^2+x(2)^2+x(3)^2)))) x(2)*x(3)/((x(1)^2+x(2)^2+x(3)^2)^(3/2)*sqrt(1-(x(3)^2/(x(1)^2+x(2)^2+x(3)^2)))) -(x(1)^2+x(2)^2)/((x(1)^2+x(2)^2+x(3)^2)^(3/2)*sqrt(1-(x(3)^2/(x(1)^2+x(2)^2+x(3)^2))))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures()    
    f2 = figure(2); clf; hold all; f2.Position = [750 200 600 475];
    view(-109,14); lighting phong; light('Position',[-1 0 0]);
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel("x", 'FontSize', 18, 'FontName', 'Times');
    ylabel("y", 'FontSize', 18, 'FontName', 'Times');
    zlabel("z", 'FontSize', 18, 'FontName', 'Times');
    %xlabel("x", 'FontSize', 18, 'FontName', 'Times', 'Position',[-10 44 -26]);
    %ylabel("y", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 -15 -42]);
    %zlabel("z", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 47 8]);
    set(get(gca,'ZLabel'), 'Rotation', 0);
    %xlim([-20 20])
    %xticks([-20 -10 0 10 20])
    %xticklabels({'-20','-10','0','10','20'})
    %ylim([-30 30])
    %yticks([-30 -20 -10 0 10 20 30])
    %yticklabels({'-30','-20','-10','0','10', '20', '30'})
    %zlim([-30 30])
    %zticks([-30 -20 -10 0 10 20 30])
    %zticklabels({'-30','-20','-10','0','10', '20', '30'})
    set(gca, 'FontName' , 'Times');
    legend('Location','Northwest')
end