clear all; close all; clc;

% PCR3BP_UKF_M.m
% Ben Hanson, 2024

%% Initial Condition
rng(2024);
po = readmatrix('./Data/periodic_orbits_JE.csv'); const.d = 4;  
const.LU = 668519; const.TU = 48562; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; 
const.T = po(9);
norm = 0; % 0: Plot is in real units, 1: plot is normalized
initialize_figures(rv,const,norm); 

%% Get Measurements
dt = 0.0001;
Q = diag([0; 0; 0; 0]);
%t = [const.T/nm const.T];
t = [const.T];
R = diag([2.2375e-08; 2.2375e-08; 5.2767e-09; 5.2767e-09]); 
[x, z, tx, tz] = get_measurements('x0', rv.start, ...
                                  'dt', dt, ...
                                  't', t, ...
                                  'f', @f, ...
                                  'h', @h, ...
                                  'Q', Q, ...
                                  'R', R, ...
                                  'const', const);

%% Plotting Nominal Trajectories
figure(1); 
plot(x(1,:).*const.LU,x(2,:).*const.LU,'Color','k','LineWidth',2,'DisplayName','Nominal');
scatter(x(1,1).*const.LU,x(2,1).*const.LU,50,'filled','MarkerFaceColor','k','HandleVisibility','off');
if(t(1) < const.T), scatter(z(1,:).*const.LU,z(2,:).*const.LU,50,'filled','MarkerFaceColor','r','DisplayName','Measurements'); end
drawnow;
figure(2); 
plot(x(3,:).*(const.LU/const.TU),x(4,:).*(const.LU/const.TU),'Color', 'k','LineWidth', 2, 'DisplayName','Nominal');
scatter(x(3,1).*(const.LU/const.TU),x(4,1).*(const.LU/const.TU),50,'filled','MarkerFaceColor','k','HandleVisibility','off');
if(t(1) < const.T), scatter(z(3,:).*(const.LU/const.TU),z(4,:).*(const.LU/const.TU),50,'filled','MarkerFaceColor','r','DisplayName','Measurements'); end
drawnow;

%% Unscented Kalman Filter
xest0 = rv.start + sqrt(Q)*randn(4,1); Pest0 = R;
filter = trackingUKF('State',xest0,...
                     'StateCovariance',Pest0,...
                     'StateTransitionFcn',@df,...
                     'ProcessNoise',Q,...
                     'MeasurementFcn',@h,...
                     'MeasurementNoise',R,...
                     'Alpha',1e-3,...
                     'Beta',2,...
                     'Kappa',0);

if length(t)==1
    m = size(R);
    z = {t,NaN(m(1),1)};
else
    if isempty(h)
        error("Included measurements but missing measurement model.")
    end
end

t0 = 0;
xest = xest0;
Pest = {}; Pest{1} = Pest0;
nz = length(tz);
t0 = 0; tspan = t0;

for i=1:nz
    tk1 = z{1}; tk1 = tk1(i); zk1 = z{2}; zk1 = zk1(:,i); t0 = tspan(end);
    
    % Prediction
    dtx = dt;
    while t0 < tk1
        dtx = min(dtx, tk1-t0); t0 = t0 + dtx; tspan(end+1) = t0;
        [xpred,Ppred] = predict(filter,dt,const);
        xest(:,end+1) = xpred;
        Pest{end+1} = Ppred;
    end

    % Correction
    if ~isnan(zk1)
        [xtemp,Ptemp] = correct(filter, zk1);
        xest(:,end) = xtemp;
        Pest{end} = Ptemp; 
    end
end

[sv, idx] = intersect(tspan,tz,"stable"); 
idx = [1; idx(2:end-1)-1; idx(end)];

figure(1);
plot(xest(1,:).*const.LU,xest(2,:).*const.LU,"m--",'LineWidth',1,"DisplayName","UKF");
scatter(xest(1,idx).*const.LU,xest(2,idx).*const.LU,50,'filled','MarkerFaceColor','m',"HandleVisibility","off");
drawnow;
figure(2); 
plot(xest(3,:).*(const.LU/const.TU),xest(4,:).*(const.LU/const.TU),"m--",'LineWidth',1,"DisplayName","UKF");
scatter(xest(3,idx).*(const.LU/const.TU),xest(4,idx).*(const.LU/const.TU),50,'filled','MarkerFaceColor','m',"HandleVisibility","off");
drawnow;

p.display = 1; 
p.name = "UKF, 3\sigma covariance";
p.color = 'm';
for i=1:length(idx)
    figure(1);
    P = Pest{idx(i)}; 
    plot_gaussian_ellipsoid(xest(1:2,idx(i)).*(const.LU),P(1:2,1:2).*(const.LU^2),3,p);
    figure(2); 
    plot_gaussian_ellipsoid(xest(3:4,idx(i)).*(const.LU/const.TU),P(3:4,3:4).*((const.LU/const.TU)^2),3,p);
    p.display = 0; 
end 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(x,const)
    
    r1 = ((x(1)+const.mu)^2+x(2)^2)^(1.5);
    r2 = ((x(1)-1+const.mu)^2+x(2)^2)^(1.5);
    
    ax = 2*x(4)+x(1)-(const.mu*(x(1)-1+const.mu)/r2)-((1-const.mu)*(x(1)+const.mu)/r1); 
    ay = -2*x(3)+x(2)-(const.mu*x(2)/r2)-((1-const.mu)*x(2)/r1);

    x1 = [x(3); x(4); ax; ay];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = df(x,dt,const)
    
    r1 = ((x(1)+const.mu)^2+x(2)^2)^(1.5);
    r2 = ((x(1)-1+const.mu)^2+x(2)^2)^(1.5);
    
    ax = 2*x(4)+x(1)-(const.mu*(x(1)-1+const.mu)/r2)-((1-const.mu)*(x(1)+const.mu)/r1); 
    ay = -2*x(3)+x(2)-(const.mu*x(2)/r2)-((1-const.mu)*x(2)/r1);

    x1 = [x(1)+dt*x(3); x(2)+dt*x(4); x(3)+dt*ax; x(4)+dt*ay];
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
function initialize_figures(rv,const,norm)

    f1 = figure(1); clf; hold all; f1.Position = [100 150 700 600]; legend;
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
    
    f2 = figure(2); clf; hold all; f2.Position = [800 150 700 600]; legend;
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