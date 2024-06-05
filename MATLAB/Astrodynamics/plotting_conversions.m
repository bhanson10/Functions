clear all; close all; clc; format long

%% Initial Conditions 
const.mu = 3202.738774922892; const.T = 10150.71332976753*3;
% rv = [1690.689890962806; 0; 0; 0; -0.2954850833646089; -1.352955033637386],
rv = [1682.559086494749; 0; 0; 0; -0.5671495124122091; -1.381021094488744],
rv_cov = diag([1^2; 1^2; 1^2; 0.001^2; 0.001^2; 0.001^2]),

%% Converting from RV and RV_cov to COE and COE_cov 
coe     = rv2oe(rv,const.mu,'rad','M','coe')',
coe_cov = rv2oe_cov(rv,rv_cov,const.mu,'M','coe'),

%% Creating RV/OE Datasets
n = 1000; 
x_rv = mvnrnd(rv,rv_cov,n); 
x_coe = mvnrnd(coe,coe_cov,n); 
for i=1:n
    x_rv_rec(i,:) = oe2rv(x_coe(i,:),const.mu,'rad', 'M', 'coe');
end

lbl.XString = '$x$ (km)'; lbl.YString = '$y$ (km)'; lbl.ZString = '$z$ (km)';
initialize_figures('n', 1, 'spacing', {[50 100 700 700]}, 'lgd', 1, 'lbl', lbl, 'vw', {[45,10]}, 'axs', {'equal'});
initialize_figures('n', 2:5, 'spacing', {[50 100 700 700]}, 'lbl', lbl, 'vw', {[45,10]});

figure(1);
scatter3(x_rv(:,1), x_rv(:,2), x_rv(:,3), 10, 'r', 'filled', 'DisplayName', 'RV');
scatter3(x_rv_rec(:,1), x_rv_rec(:,2), x_rv_rec(:,3), 10, 'b', 'filled', 'DisplayName', 'OE');
drawnow; 

%% Propagating OEs
options = odeset('MaxStep', 5, 'InitialStep', 0.1, 'RelTol', 1e-6);
for i = 1:n
    x0 = x_coe(i,:); 
    [ti, xi] = ode87(@(t, x) f(t, x, const), [0,const.T], x0, options);
    nt = length(ti); 
    x_rv_T25(i,:) = oe2rv(xi(round(nt/4),:), const.mu,'rad', 'M', 'coe');
    x_rv_T50(i,:) = oe2rv(xi(round(nt/2),:), const.mu,'rad', 'M', 'coe');
    x_rv_T75(i,:) = oe2rv(xi(round(3*nt/4),:), const.mu,'rad', 'M', 'coe');
    x_rv_T(i,:) = oe2rv(xi(end,:), const.mu,'rad', 'M', 'coe');
end

figure(2);
scatter3(x_rv_T25(:,1), x_rv_T25(:,2), x_rv_T25(:,3), 10, 'b', 'filled', 'HandleVisibility', 'off');
drawnow; 
figure(3);
scatter3(x_rv_T50(:,1), x_rv_T50(:,2), x_rv_T50(:,3), 10, 'b', 'filled', 'HandleVisibility', 'off');
drawnow; 
figure(4);
scatter3(x_rv_T75(:,1), x_rv_T75(:,2), x_rv_T75(:,3), 10, 'b', 'filled', 'HandleVisibility', 'off');
drawnow; 
figure(5);
scatter3(x_rv_T(:,1), x_rv_T(:,2), x_rv_T(:,3), 10, 'b', 'filled', 'HandleVisibility', 'off');
drawnow; 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(t, x,const)
    x1 = [0; 0; 0; 0; 0; sqrt(const.mu/x(1)^3)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%