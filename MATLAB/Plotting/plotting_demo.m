clear all; clc; close all; 

mu = [0; 0; 0]; S = [1 0 0; 0 1 0; 0 0 1]; 

n = 100000; 
x = mvnrnd(mu,S,n);

for i=1:n, P(i) = exp(-0.5*(x(i,:)-mu')*inv(S)*(x(i,:)-mu')'); end, P = P./sum(P); 
colors = jet(100);
scaled_P = (P - min(P))./(max(P)-min(P));
colors = colors(round(scaled_P.*99)+1,:);

figure(1); hold on; view(45,20); axis equal; title('PF: Point Mass');  
scatter3(x(:,1), x(:,2), x(:,3), 10, colors, 'filled'); 

figure(2); hold on; view(45,20); axis equal; title('PF: Isosurface'); 
plot_nongaussian_surface(x,P,[0.2, 0.74, 0.971], 0.1);

figure(3); hold on; view(45,20); axis equal; title('Ellipsoids');  
plot_gaussian_ellipsoid(mu,S,[1, 2, 3]);
