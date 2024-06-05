clear all; close all; clc;

% PCR3BP_EKF.m
% Benjamin Hanson, 2024
%
% In this example, an EKF performs state estimation on a space object with
% initial conditions and equations of motion representative of a L3
% Lyapunov orbit of the Jupiter-Europa System

po = readmatrix('./Data/periodic_orbits_JE.csv'); const.d = 4;  
const.LU = 668519; const.TU = 48562; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; 
const.T = po(9);
norm = 0; % 0: Plot is in real units, 1: plot is normalized
initialize_figures(rv,const,norm); 

rng(2024);                                        % Random Number Generator
dt = 0.0001;                                 % Measurement cadence
Q = diag([0; 0; 0; 0]);                           % Process noise matrix
R = diag([2.25E-6; 2.25E-6; 5.2E-5; 5.2E-5]);     % Measurement noise matrix. Units are LU and LU/TU.

% True States and Measurements
[x, z, tx, tz] = get_truth_measurements('x0', rv.start, ...
                                        'dt', dt(1), ...
                                        't', [const.T], ...  % a.k.a no measurements
                                        'f', @f, ...
                                        'h', @h, ...
                                        'Q', Q, ...
                                        'R', R, ...
                                        'const', const, ...
                                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation

% Extended Kalman Filter
xest0 = rv.start + sqrt(R)*randn(4,1); Pest0 = R;
[xest,Pest,txest] = EKF('xest0', xest0, ...
                        'Pest0', Pest0, ...
                        'dt', dt, ...
                        'z', {const.T}, ...  % a.k.a no measurements
                        'Q', Q, ...
                        'R', R, ...
                        'f', @f, ...
                        'F', @F, ...
                        'const', const, ...
                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation/EKF

figure(1)
subplot(4,1,1)
plot(tz,z(1,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("x (LU)")
subplot(4,1,2)
plot(tz,z(2,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("y (LU)")
subplot(4,1,3)
plot(tz,z(3,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("v_x (LU/TU)")
subplot(4,1,4)
plot(tz,z(4,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("v_y (LU/TU)")

% Plotting MC Trajectories
n = 200; 
YX = []; YY = []; YVX = []; YVY = []; YY0 = zeros(4,n);
for i = 1:n
    Y0 = rv.start + sqrt(R)*randn(length(R),1); YY0(:,i) = Y0; tspan = [0 const.T]; 
    options = odeset('RelTol', 1e-13); % Setting a tolerance
    [t, Y] = ode45(@(t, Y) f(Y,const), tspan, Y0, options);
    
    if(norm)
        figure(2); 
        plot(Y(:,1),Y(:,2),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off');
        YX(end+1) = Y(end,1); YY(end+1) = Y(end,2);
        drawnow; 
        figure(3); 
        plot(Y(:,3),Y(:,4),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off');
        YVX(end+1) = Y(end,3); YVY(end+1) = Y(end,4);
        drawnow; 
    else
        figure(2); 
        plot(Y(:,1).*const.LU,Y(:,2).*const.LU,'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off');
        YX(end+1) = Y(end,1).*const.LU; YY(end+1) = Y(end,2).*const.LU;
        drawnow; 
        figure(3);
        plot(Y(:,3).*(const.LU/const.TU),Y(:,4).*(const.LU/const.TU),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off');
        YVX(end+1) = Y(end,3).*(const.LU/const.TU); YVY(end+1) = Y(end,4).*(const.LU/const.TU);
        drawnow;
    end
end
figure(2); 
scatter(YY0(1,:),YY0(2,:),50,'k','filled', 'HandleVisibility','off');
scatter(YX,YY,20,'k','filled', 'HandleVisibility','off');
figure(3); 
scatter(YY0(3,:),YY0(4,:),50,'k','filled', 'HandleVisibility','off');
scatter(YVX,YVY,20,'k','filled', 'HandleVisibility','off');

if(norm)
    figure(2); 
    cad = round(linspace(1,length(txest),10));
    plot(x(1,:),x(2,:),"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(1,:),xest(2,:),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(1,cad),x(2,cad),'r', 'filled','HandleVisibility','off');
    scatter(xest(1,cad),xest(2,cad),'b', 'filled','HandleVisibility','off');
    drawnow;
    figure(3); 
    cad = round(linspace(1,length(tx),10));
    plot(x(3,:),x(4,:),"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(3,:),xest(4,:),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(3,cad),x(4,cad),'r', 'filled','HandleVisibility','off');
    scatter(xest(3,cad),xest(4,cad),'b', 'filled','HandleVisibility','off');
    drawnow;
else
    figure(2);
    cad = round(linspace(1,length(tx),10));
    plot(x(1,:).*const.LU,x(2,:).*const.LU,"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(1,:).*const.LU,xest(2,:).*const.LU,"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(1,cad).*const.LU,x(2,cad).*const.LU,'r', 'filled','HandleVisibility','off');
    scatter(xest(1,cad).*const.LU,xest(2,cad).*const.LU,'b', 'filled','HandleVisibility','off');
    drawnow;
    figure(3); 
    cad = round(linspace(1,length(tx),10));
    plot(x(3,:).*(const.LU/const.TU),x(4,:).*(const.LU/const.TU),"r-",'LineWidth',1,"DisplayName","True Trajectory");
    plot(xest(3,:).*(const.LU/const.TU),xest(4,:).*(const.LU/const.TU),"b--",'LineWidth',1,"DisplayName","Estimated Trajectory");
    scatter(x(3,cad).*(const.LU/const.TU),x(4,cad).*(const.LU/const.TU),'r', 'filled','HandleVisibility','off');
    scatter(xest(3,cad).*(const.LU/const.TU),xest(4,cad).*(const.LU/const.TU),'b', 'filled','HandleVisibility','off');
    drawnow;
end

p.display = 1; 
p.name = "3\sigma covariance";
for i=round(linspace(1,length(tx),10))
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
function x1 = f(X,const)     
    x = X(1); y = X(2); vx = X(3); vy = X(4); 

    r1 = ((x+const.mu)^2+y^2)^(1.5);
    r2 = ((x-1+const.mu)^2+y^2)^(1.5);
    
    ax = 2*vy+x-(const.mu*(x-1+const.mu)/r2)-((1-const.mu)*(x+const.mu)/r1); 
    ay = -2*vx+y-(const.mu*y/r2)-((1-const.mu)*y/r1);

    x1 = [vx; vy; ax; ay]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x)
    z = [x(1); x(2); x(3); x(4)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = F(X,const)
    x = X(1); y = X(2);

    
    Fx =   [0 0 1 0;
            0 0 0 1;
            ((const.mu - 1)/((const.mu + x)^2 + y^2)^(3/2) - const.mu/((const.mu + x - 1)^2 + y^2)^(3/2) + (3*const.mu*(2*const.mu + 2*x - 2)*(const.mu + x - 1))/(2*((const.mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*const.mu + 2*x)*(const.mu + x)*(const.mu - 1))/(2*((const.mu + x)^2 + y^2)^(5/2)) + 1) -((3*y*(const.mu + x)*(const.mu - 1))/((const.mu + x)^2 + y^2)^(5/2) - (3*const.mu*y*(const.mu + x - 1))/((const.mu + x - 1)^2 + y^2)^(5/2)) 0 2;
            ((3*const.mu*y*(2*const.mu + 2*x - 2))/(2*((const.mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*const.mu + 2*x)*(const.mu - 1))/(2*((const.mu + x)^2 + y^2)^(5/2))) ((const.mu - 1)/((const.mu + x)^2 + y^2)^(3/2) - const.mu/((const.mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(const.mu - 1))/((const.mu + x)^2 + y^2)^(5/2) + (3*const.mu*y^2)/((const.mu + x - 1)^2 + y^2)^(5/2) + 1) -2 0];
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
function [value, isterminal, direction] = xcross_event(~, X)
% X - Crossing Event Function
value = X(2);
isterminal=  0; %   1 = end integration
                %   0 = continue integration
direction = -1; %   1 = crossing with ydot > 0
                %  -1 = crossing with ydot < 0
                %   0 = doesn't matter (includes all)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%