clear all; clc; close all; 

load('3D_scatter_X.mat');
load('3D_scatter_P.mat');
colors = jet(100);
scaled_P = (P - min(P))./(max(P)-min(P));
colors = colors(round(scaled_P.*99)+1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lbl.String = {"x (km)", "z (km)"}; lbl.FontSize=16; lbl.FontName = 'Times';
lmts.XLimit = [min(X(:,1)),max(X(:,1))]; lmts.YLimit = [min(X(:,3)),max(X(:,3))];
initialize_figures('n', 1:3, 'spacing', {[50 100 700 700],[750 100 700 700],[750 100 700 700]}, 'lbl', lbl, 'lmts', lmts);

figure(1); 
scatter(X(:,1),X(:,3), 5, colors, 'filled');

figure(2);   
p.color = 'blue'; p.means = 1; p.type = 'smooth';
[~,shps] = plot_nongaussian_surface(X(:,[1,3]),P,[0.68,0.957,0.997],0.2,p);

figure(3); 
scatter(X(:,1),X(:,3), 5, colors, 'filled');
[~,shps] = plot_nongaussian_surface(X(:,[1,3]),P,[0.68,0.957,0.997],0.2,p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lbl.String = {"x (km)", "y (km)", "z (km)"}; lbl.FontSize=16; lbl.FontName = 'Times';
lmts.XLimit = [min(X(:,1)),max(X(:,1))]; lmts.YLimit = [min(X(:,2)),max(X(:,2))]; lmts.ZLimit = [min(X(:,3)),max(X(:,3))];
vw.az = 45; vw.el = 30; 
initialize_figures('n', 1:3, 'spacing', {[50 100 700 700],[750 100 700 700],[750 100 700 700]}, 'lbl', lbl, 'lmts', lmts, 'vw', vw);

figure(1); 
scatter3(X(:,1),X(:,2),X(:,3), 5, colors, 'filled');

figure(2);
p.color = 'blue'; p.means = 1; p.type = 'smooth';
[~,shps] = plot_nongaussian_surface(X,P,[0.68,0.957,0.997],0.2,p);

figure(3);
scatter3(X(:,1),X(:,2),X(:,3), 5, colors, 'filled');
[~,shps] = plot_nongaussian_surface(X,P,[0.68,0.957,0.997],0.2,p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}