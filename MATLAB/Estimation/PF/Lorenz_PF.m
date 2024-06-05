clear all; close all; clc;

% Lorenz_PF.m
% Benjamin Hanson, 2024
%
% Applying the PF to the Lorenz system

% Initial conditions
initialize_figures(); 
rng(2024); % Random Number Generator
const.sigma=4; const.b=1; const.r=48; % Lorenz Parameters
T   = 1; % Period of propagation
dt  = 0.02; % Time step or [True time step, Measurement time step]
x0  = [-11.5; -10; 9.5]; % [x;y;z]
Q   = diag([.01; .01; .01;]); % Process noise matrix
R  = diag([1; 0.1; 0.1]); % Measurement noise matrix, spherical model 
%R   = [1]; % Measurement noise matrix, z-value only

% True States and Measurements
[x, z, tx, tz] = get_truth_measurements('x0', x0, ...
                                        'dt', dt, ...
                                        't', [dt,T], ...
                                        'f', @f, ...
                                        'h', @h, ...
                                        'Q', Q, ...
                                        'R', R, ...
                                        'const', const, ...
                                        'method', 'RK4', ...
                                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation

% Particle Filter
num = 10000;
[xk,wk,xest,txest] = PF('xest0', x0, ...
                        'Pest0', diag([1,1,1]), ...
                        'dt', dt, ...
                        'z', {tz, z}, ... 
                        'Q', Q, ...
                        'R', R, ...
                        'f', @f, ...
                        'h', @h, ...
                        'l', @l, ...
                        'num', num, ...
                        'const', const,...
                        'filter', 'SIS', ...
                        'time', 'CT'); %$.../Utilities/Functions/MATLAB/Estimation/PF

% Plotting Background behavior
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) f(Y,const), [0,50], x0, options);
figure(2); 
plot3(Y(:,1),Y(:,2),Y(:,3),'Color', '[0 1 0 0.5]','linewidth',0.3, 'HandleVisibility','off');
drawnow;

figure(1)
subplot(3,1,1)
plot(tz,z(1,:),'LineWidth',1)
xlabel("time (TU)")
ylabel("range (LU)")
subplot(3,1,2)
plot(tz,z(2,:)*180/pi,'LineWidth',1)
xlabel("time (TU)")
ylabel("\theta (deg)")
subplot(3,1,3)
plot(tz,z(3,:)*180/pi,'LineWidth',1)
xlabel("time (TU)")
ylabel("\phi (deg)")

% figure(1)
% subplot(1,1,1)
% plot(tz,z(1,:),'LineWidth',1)
% xlabel("time (s)")
% ylabel("z (LU)")

figure(2); hold on;
cad = round(linspace(1,length(tx),6));
plot3(x(1,:),x(2,:),x(3,:),"k-",'LineWidth',1,"DisplayName","True Trajectory");
scatter3(x(1,cad),x(2,cad),x(3,cad),'k', 'filled','HandleVisibility','off');
for i=cad
    x = []; w = []; colors = jet(100);
    for j=1:num
        xj = xk{j}; wj = wk{j};
        x(:,end+1) = xj(:,i); w(end+1) = wj(i);
    end
    if all(w == w(1))
        % If all weights are the same, assign a single color
        colors = repmat([0 0 0.52], num, 1);
    else
        % Define colors based on weights
        scaled_w = (w - min(w)) / (max(w) - min(w)); % Normalize weights to [0,1]
        colors = colors(round(scaled_w * 99) + 1,:); % Map weights to colors
    end
    scatter3(x(1,:),x(2,:),x(3,:), 5, colors, 'filled','MarkerFaceAlpha', 0.2, 'HandleVisibility','off');
    %scatter3(xest(1,i),xest(2,i),xest(3,i), 50, 'd', 'm', 'filled','HandleVisibility','off');
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=f(x,const)                          
    f=[const.sigma*(x(2)-x(1));  -x(2)-x(1)*x(3);  -const.b*x(3)+x(1)*x(2)-const.b*const.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=l(x,mu,Sigma,const)
    output = exp(-(0.5*(x-mu)'*inv(Sigma)*(x-mu)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x) % spherical model
    z = [sqrt(x(1)^2+x(2)^2+x(3)^2); atan2(x(2),x(1)); acos(x(3)/sqrt(x(1)^2+x(2)^2+x(3)^2))]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h2(x) % z-value model
    z = [x(3)]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = F(x,const)
    Fx = [-const.sigma const.sigma 0       ;         
          -x(3)        -1          -x(1)   ;
          x(2)         x(1)        -const.b];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%