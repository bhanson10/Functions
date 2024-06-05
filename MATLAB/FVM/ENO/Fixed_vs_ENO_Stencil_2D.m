clear all; close all; clc; 

%% Defining x_i and v_i
v = @(x,y,const) double(max(0,1-sqrt((x-const.c_c(1)).^2+(y-const.c_c(2)).^2)*(1/const.r_c)))+double((abs(x-const.c_wtc(1))<=const.r_wtc)&(abs(y-const.c_wtc(2))<=const.r_wtc));
% v=@(x,y,const) double((x+y)<=const.mu)+double((x+y)>const.mu)/2; % 2. Diagonal Step Function

const.mu=1; const.sig=0.5; const.c_wtc=[0.4,0]; const.r_wtc=0.3; const.c_c=[-0.4,0]; const.r_c=0.3;
xl=[-1, 1]; yl = [-1, 1];
Nx=201; dx=(xl(2)-xl(1))/(Nx-1); xbar=xl(1):dx:xl(2);  
Ny=201; dy=(yl(2)-yl(1))/(Ny-1); ybar=yl(1):dy:yl(2);  
vbar=NaN(Nx,Ny); 
[Xbar,Ybar]=meshgrid(xbar,ybar);

options = {'Method','iterated','AbsTol',1e-6,'RelTol',1e-4};
for i=1:Nx
    for j=1:Ny
        % xim = xbar(i)-(dx/2); xip = xbar(i)+(dx/2);
        % yim = ybar(j)-(dy/2); yip = ybar(j)+(dy/2);
        % vbar(i,j) = (1/(dx*dy))*integral2(@(x,y) v(x,y,const), xim, xip, yim, yip, options{:}); % Cell Averages
        vbar(i,j)=v(xbar(i),ybar(j),const);
    end
end

%% Plotting
lbl.String = {'$x$', '$y$'}; lbl.FontSize = 24; lbl.Interpreter = 'latex'; 
initialize_figures('n', 1, 'margin', [800 200]', 'lbl', lbl, 'vw', {[45, 30]}, 'lgd', 1)
mesh(Xbar,Ybar,vbar,'EdgeColor','r','FaceColor','r','DisplayName','$\bar{v}_{ij}$'); alpha(.2);
