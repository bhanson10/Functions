clear all; close all; clc;

% PCR3BP_EKF_M.m
% Ben Hanson
% 1/16/2024 
%
% In this example, an EKF performs state estimation on a space object with
% initial conditions and equations of motion representative of a L3
% Lyapunov orbit of the Jupiter-Europa System using a built-in MATLAB
% function

po = readmatrix('./Data/periodic_orbits_JE.csv'); const.d = 4;  
const.LU = 668519; const.TU = 48562; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; 
const.T = po(9);
norm = 1; % 0: Plot is in real units, 1: plot is normalized
initialize_figures(rv,const,norm); 

rng(2024);                                        % Random Number Generator
dt = [0.0001,10];                                 % Measurement cadence
Q = diag([0; 0; 0; 0]);                           % Process noise matrix
R = diag([2.25E-6; 2.25E-6; 5.2E-5; 5.2E-5]);     % Measurement noise matrix. Units are LU and LU/TU.
xest0 = rv.start + sqrt(Q)*randn(4,1); Pest0 = R; % Initial Estimation

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
    tspan_x = 0:dt:const.T;
    tspan_z = 0:dt:const.T;
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
    tspan_x = 0:dt(1):const.T;
    tspan_z = 0:dt(2):const.T;
    dt_ratio = dt(2)/dt(1);
    dt = dt(1); 
else
    error('Input time step must be of length 1 or 2.');
end

d1  = length(Q); d2 = length(R); 

n1 = length(tspan_x); n2 = length(tspan_z);
x  = NaN(d1,n1); x(:,1) = rv.start;
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

figure(1)
subplot(4,1,1)
plot(tspan_z,z(1,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("x (LU)")
subplot(4,1,2)
plot(tspan_z,z(2,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("y (LU)")
subplot(4,1,3)
plot(tspan_z,z(3,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("v_x (LU/TU)")
subplot(4,1,4)
plot(tspan_z,z(4,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("v_y (LU/TU)")

%{
% Plotting nominal trajectories
Y0 = rv.start; tspan = [0 const.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) PCR3BP(Y,const), tspan, Y0, options);

if(norm)
    figure(2); 
    plot(Y(:,1),Y(:,2),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    figure(3); 
    plot(Y(:,3),Y(:,4),'k-','linewidth',1.5,'DisplayName','Nominal');
    drawnow; 
else
    figure(2); 
    plot(Y(:,1).*const.LU,Y(:,2).*const.LU,'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    figure(3);
    plot(Y(:,3).*(const.LU/const.TU),Y(:,4).*(const.LU/const.TU),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow;
end
%}

if(norm)
    figure(2); 
    cad = round(linspace(1,length(tspan_x),10));
    plot(x(1,:),x(2,:),"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(1,:),xest(2,:),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(1,cad),x(2,cad),'r', 'filled','HandleVisibility','off');
    scatter(xest(1,cad),xest(2,cad),'b', 'filled','HandleVisibility','off');
    drawnow;
    figure(3); 
    cad = round(linspace(1,length(tspan_x),10));
    plot(x(3,:),x(4,:),"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(3,:),xest(4,:),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(3,cad),x(4,cad),'r', 'filled','HandleVisibility','off');
    scatter(xest(3,cad),xest(4,cad),'b', 'filled','HandleVisibility','off');
    drawnow;
else
    figure(2);
    cad = round(linspace(1,length(tspan_x),10));
    plot(x(1,:)*const.LU,x(2,:)*const.LU,"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(1,:).*const.LU,xest(2,:).*const.LU,"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(1,cad).*const.LU,x(2,cad).*const.LU,'r', 'filled','HandleVisibility','off');
    scatter(xest(1,cad).*const.LU,xest(2,cad).*const.LU,'b', 'filled','HandleVisibility','off');
    drawnow;
    figure(3); 
    cad = round(linspace(1,length(tspan_x),10));
    plot(x(3,:),x(4,:),"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(3,:).*(const.LU/const.TU),xest(4,:).*(const.LU/const.TU),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(3,cad).*(const.LU/const.TU),x(4,cad).*(const.LU/const.TU),'r', 'filled','HandleVisibility','off');
    scatter(xest(3,cad).*(const.LU/const.TU),xest(4,cad).*(const.LU/const.TU),'b', 'filled','HandleVisibility','off');
    drawnow;
end

p.display = 1; 
p.name = "3\sigma covariance";
for i=round(linspace(1,length(tspan_x),10))
    if(norm)
        figure(2);
        P = Pest{i}; 
        plot_gaussian_ellipsoid(xest(1:2,i),P(1:2,1:2),3,p);
        figure(3); 
        plot_gaussian_ellipsoid(xest(3:4,i),P(3:4,3:4),3,p);
    else
        figure(2);
        P = Pest{i}; 
        plot_gaussian_ellipsoid(xest(1:2,i).*(const.LU),P(1:2,1:2).*(const.LU^2),3,p);
        figure(3); 
        plot_gaussian_ellipsoid(xest(3:4,i).*(const.LU/const.TU),P(3:4,3:4).*((const.LU/const.TU)^2),3,p);
    end
    p.display = 0; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=PCR3BP(X,const)     
    x = X(1); y = X(2); vx = X(3); vy = X(4); 

    r1 = ((x+const.mu)^2+y^2)^(1.5);
    r2 = ((x-1+const.mu)^2+y^2)^(1.5);
    
    ax = 2*vy+x-(const.mu*(x-1+const.mu)/r2)-((1-const.mu)*(x+const.mu)/r1); 
    ay = -2*vx+y-(const.mu*y/r2)-((1-const.mu)*y/r1);

    f = [vx; vy; ax; ay]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(x,dt,const)
    
    r1 = ((x(1)+const.mu)^2+x(2)^2)^(1.5);
    r2 = ((x(1)-1+const.mu)^2+x(2)^2)^(1.5);
    
    ax = 2*x(4)+x(1)-(const.mu*(x(1)-1+const.mu)/r2)-((1-const.mu)*(x(1)+const.mu)/r1); 
    ay = -2*x(3)+x(2)-(const.mu*x(2)/r2)-((1-const.mu)*x(2)/r1);

    x1 = [x(1) + x(3)*dt; x(2) + x(4)*dt; x(3) + ax*dt; x(4) + ay*dt];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x)
    z = [x(1); x(2); x(3); x(4)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = F(X,dt,const)
    x = X(1); y = X(2);

    
    Fx =   [1 0 dt 0;
            0 1 0 dt;
            dt*((const.mu - 1)/((const.mu + x)^2 + y^2)^(3/2) - const.mu/((const.mu + x - 1)^2 + y^2)^(3/2) + (3*const.mu*(2*const.mu + 2*x - 2)*(const.mu + x - 1))/(2*((const.mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*const.mu + 2*x)*(const.mu + x)*(const.mu - 1))/(2*((const.mu + x)^2 + y^2)^(5/2)) + 1) -dt*((3*y*(const.mu + x)*(const.mu - 1))/((const.mu + x)^2 + y^2)^(5/2) - (3*const.mu*y*(const.mu + x - 1))/((const.mu + x - 1)^2 + y^2)^(5/2)) 1 2*dt;
            dt*((3*const.mu*y*(2*const.mu + 2*x - 2))/(2*((const.mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*const.mu + 2*x)*(const.mu - 1))/(2*((const.mu + x)^2 + y^2)^(5/2))) dt*((const.mu - 1)/((const.mu + x)^2 + y^2)^(3/2) - const.mu/((const.mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(const.mu - 1))/((const.mu + x)^2 + y^2)^(5/2) + (3*const.mu*y^2)/((const.mu + x - 1)^2 + y^2)^(5/2) + 1) -2*dt 1];
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hx = H(x)
    Hx = [1 0 0 0;         
          0 1 0 0;
          0 0 1 0;
          0 0 0 1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = plotErrorEllipse(P, p)
    s = -2 * log(1 - p);
    [V, D] = eig(P * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures(rv,const,norm)

    f2 = figure(2); clf; hold all; f2.Position = [100 150 700 600]; legend;
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlabel("x (LU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (LU)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-1.75, 1.25])
        ylim([-1, 1])
        scatter(-const.mu,0,75,'filled','MarkerFaceColor','m','HandleVisibility','off');
        text(-const.mu,-0.08, 'Jupiter', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter(1-const.mu,0,30,'filled','MarkerFaceColor','#808080','HandleVisibility','off');
        text(1-const.mu,-0.08, 'Europa', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter(-1.00001053,0,50,'d','k','HandleVisibility','off')
        text(-0.9,-0.08, 'L_3','HorizontalAlignment','center','FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        %scatter(rv.start(1),rv.start(2),25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    else
        xlabel("x (km)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (km)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-1.75*const.LU, 1.25*const.LU])
        ylim([-1*const.LU, 1*const.LU])
        scatter(-const.mu*const.LU,0,75,'filled','MarkerFaceColor','m','HandleVisibility','off');
        text(-const.mu*const.LU,-0.08*const.LU, 'Jupiter', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter((1-const.mu)*const.LU,0,30,'filled','MarkerFaceColor','#808080','HandleVisibility','off');
        text((1-const.mu)*const.LU,-0.08*const.LU, 'Europa', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter(-1.00001053*const.LU,0,50,'d','k','HandleVisibility','off')
        text(-0.9*const.LU,-0.08*const.LU, 'L_3','HorizontalAlignment','center','FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        %scatter(rv.start(1)*const.LU,rv.start(2)*const.LU,25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    end
    drawnow; 
    
    f3 = figure(3); clf; hold all; f3.Position = [800 150 700 600]; legend;
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlabel("v_x (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.85, 0.85])
        ylim([-1.2, 1])
        %scatter(rv.start(3),rv.start(4),25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    else
        xlabel("v_x (km/s)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (km/s)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.85*(const.LU/const.TU), 0.85*(const.LU/const.TU)])
        ylim([-1.2*(const.LU/const.TU), 1*(const.LU/const.TU)])
        %scatter(rv.start(3)*(const.LU/const.TU),rv.start(4)*(const.LU/const.TU),25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    end
    drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%