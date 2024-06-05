clear all; close all; clc;

% Lorenz_EKF_M.m
% Ben Hanson, 2024
%
% Applying the EKF to the Lorenz system, using built-in MATLAB function

% Initial conditions
initialize_figures(); 
rng(3); % Random Number Generator
const.sigma=4; const.b=1; const.r=48; % Lorenz Parameters
T   = 1; % Period of propagation
dt  = [0.02,2]; % Time step or [True time step, Measurement time step]
x0  = [-11.5; -10; 9.5]; % [x;y;z]
Q   = diag([.01; .01; .01;]); % Process noise matrix
R   = diag([1; 0.1; 0.1]); % Measurement noise matrix
xest0 = x0+sqrt(R(1,1))*randn(3,1); Pest0 = diag([1,1,1]); % Initial Estimate

% Extended Kalman Filter
filter = trackingEKF("State", xest0,...
                     "StateCovariance",Pest0,...
                     "StateTransitionFcn",@f,...
                     "StateTransitionJacobianFcn",@F,...
                     "ProcessNoise",Q,...
                     "MeasurementFcn",@h,...
                     "MeasurementJacobianFcn",@H,...
                     "MeasurementNoise",R);

if(length(dt)==1)
    tspan_x = 0:dt:T;
    tspan_z = 0:dt:T;
    dt_ratio = 1;
elseif(length(dt)==2)
    if(dt(2) < dt(1))
        disp('Measurment time step must be larger than true time step. Setting measurement timestep equal to true time step.')
        dt(2) = dt(1);   
    else
        if(mod(dt(2)/dt(1),1)~=0)
            disp("Measurement time step is not a multiple of true time step. Rounding to nearest multiple.")
            dt(2) = round(dt(2)/dt(1))*dt(1); 
        end
    end
    tspan_x = 0:dt(1):T;
    tspan_z = 0:dt(2):T;
    dt_ratio = dt(2)/dt(1);
    dt = dt(1); 
else
    error('Input time step must be of length 1 or 2.');
end

d1  = length(Q); d2 = length(R); 

n1 = length(tspan_x); n2 = length(tspan_z);
x  = NaN(d1,n1); x(:,1) = x0;
xest = NaN(d1,n1); xest(:,1) = xest0;
Pest = {}; Pest{1} = Pest0;
z  = NaN(d2,n2); count = 1; 

for i = 2:n1
    
    % Initialization
    w = sqrt(Q)*randn(d1,1); 
    x(:,i) = f(x(:,i-1),dt,const) + w; % True state

    % Prediction
    [xpred,Ppred] = predict(filter, dt, const);

    % Correction
    if(mod(i-1,dt_ratio)==0)
        v = sqrt(R)*randn(d2,1);
        z(:,count) = h(x(:,i)) + v; % Measurement
        [xtemp,Ptemp] = correct(filter,z(:,count));  
        xest(:,i) = xtemp;
        Pest{i} = Ptemp;
        count = count + 1; 
    else
        xest(:,i) = xpred;
        Pest{i} = Ppred; 
    end
end

% Plotting Background behavior
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,const), [0,50], x0, options);
figure(2); 
plot3(Y(:,1),Y(:,2),Y(:,3),'Color', '[0 1 0 0.5]','linewidth',0.3, 'HandleVisibility','off');
drawnow;

figure(1)
subplot(3,1,1)
plot(tspan_z,z(1,:),'LineWidth',1)
xlabel("time (s)")
ylabel("range (m)")
subplot(3,1,2)
plot(tspan_z,z(2,:)*180/pi,'LineWidth',1)
xlabel("time (s)")
ylabel("\theta (deg)")
subplot(3,1,3)
plot(tspan_z,z(3,:)*180/pi,'LineWidth',1)
xlabel("time (s)")
ylabel("\phi (deg)")

figure(2); hold on;
cad = round(linspace(1,length(tspan_x),6));
plot3(x(1,:),x(2,:),x(3,:),"m-",'LineWidth',1,"DisplayName","True Trajectory");
plot3(xest(1,:),xest(2,:),xest(3,:),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
scatter3(x(1,cad),x(2,cad),x(3,cad),'m', 'filled','HandleVisibility','off');
scatter3(xest(1,cad),xest(2,cad),xest(3,cad),'b', 'filled','HandleVisibility','off');

for i=round(linspace(1,length(tspan_x),6))
    plot_gaussian_ellipsoid(xest(1:3,i),Pest{i},3,i);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=Lorenz3D(y,const)                          
    f=[const.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -const.b*y(3)+y(1)*y(2)-const.b*const.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(x,dt,const)
    x1 = [x(1)+(const.sigma*(x(2)-x(1)))*dt;  x(2)+(-x(2)-x(1)*x(3))*dt;  x(3)+(-const.b*x(3)+x(1)*x(2)-const.b*const.r)*dt]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x)
    z = [sqrt(x(1)^2+x(2)^2+x(3)^2); atan(x(2)/x(1)); acos(x(3)/sqrt(x(1)^2+x(2)^2+x(3)^2))]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = F(x,dt,const)
    Fx = [1-const.sigma*dt const.sigma*dt 0;         
          -x(3)*dt         1-dt           -x(1)*dt;
          x(2)*dt          x(1)*dt        1-const.b*dt];
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = H(x)
    Hx = [x(1)/sqrt(x(1)^2+x(2)^2+x(3)^2) x(2)/sqrt(x(1)^2+x(2)^2+x(3)^2) x(3)/sqrt(x(1)^2+x(2)^2+x(3)^2);
          -x(2)/(x(1)^2+x(2)^2)           x(1)/(x(1)^2+x(2)^2)            0;
          x(1)*x(3)/((x(1)^2+x(2)^2+x(3)^2)^(3/2)*sqrt(1-(x(3)^2/(x(1)^2+x(2)^2+x(3)^2)))) x(2)*x(3)/((x(1)^2+x(2)^2+x(3)^2)^(3/2)*sqrt(1-(x(3)^2/(x(1)^2+x(2)^2+x(3)^2)))) -(x(1)^2+x(2)^2)/((x(1)^2+x(2)^2+x(3)^2)^(3/2)*sqrt(1-(x(3)^2/(x(1)^2+x(2)^2+x(3)^2))))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = plotErrorEllipse(P, p)
    s = -2 * log(1 - p);
    [V, D] = eig(P * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures()    
    f2 = figure(2); %clf; hold all; f2.Position = [750 200 600 475];
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