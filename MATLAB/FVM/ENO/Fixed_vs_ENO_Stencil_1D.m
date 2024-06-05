clear all; close all; clc; 

%% Defining x_i and v_i 
v=@(x, const) (x<=const.mu)-(x>const.mu); % 1. Step Function
% v=@(x,const) exp(-((x-const.mu)/const.sig).^2); % 2. Gaussian
% v=@(x,const) (x<=const.mu).*exp(-((x-const.mu)/const.sig).^2)-(x>const.mu).*exp(-((x-const.mu)/const.sig).^2); % 3. Gaussian + Step

const.mu=0; const.sig=0.5;
xl=[-1, 1]; N=10; dx=(xl(2)-xl(1))/(N-1); xbar=xl(1):dx:xl(2);  
vbar=NaN(N,1); 
for i=1:N
    xim = xbar(i)-(dx/2); xip = xbar(i)+(dx/2);
    vbar(i) = (1/dx)*integral(@(x) v(x, const), xim, xip); % Cell Averages
end

%% Plotting
lbl.String = {'$x$', '$v$'}; lbl.FontSize = 24; lbl.Interpreter = 'latex';  
initialize_figures('n', 1, 'margin', [800 200]', 'lgd', 1, 'lbl', lbl);

figure(1); 
plot(linspace(xl(1),xl(2),100000),v(linspace(xl(1),xl(2),100000),const),'r-','LineWidth',1,'DisplayName','$v(x)$');
bar(xbar,vbar,1,'FaceAlpha',0.2,'DisplayName','$\bar{v}_i$');

%% Fixed vs. Adaptive stencil
k = 4; V = ND1D(xbar, vbar, dx, k); xp = xbar+dx/2; ip = 1:N;

% Generating fixed stencil vector
if mod(k,2)==0, r = ones(1,N).*((k/2)-1); else, r = ones(1,N).*((k-1)/2); end
r(r>ip-1)=ip(r>ip-1)-1; r(r<ip+k-1-N)=k-(N-ip(r<ip+k-1-N)+1);

[vf, Sf] = F1DNR(xbar,V,dx,xp,ip,k,r); % Fixed Stencil
[va, Sa] = ENO1DNR(xbar,V,dx,xp,ip,k); % ENO Stencil

plot(xp, vf, 'go-','LineWidth',1,'DisplayName','Fixed');
plot(xp, va, 'bo-','LineWidth',1,'DisplayName','Adaptive');

f2 = figure(2); clf; hold on; f2.Position = [50 100 1000 400];  
tcl = tiledlayout(1,2);

nexttile(tcl); hold all; legend('Location', 'southeast', 'interpreter','latex');
set(gca, 'FontName' , 'Times','FontSize',14); 
xlabel('$x$', 'Interpreter','latex', 'FontSize',24);
ylabel('$I_i$', 'Interpreter','latex', 'FontSize',24);
ylim([0.5,N+0.5]);
yticks(1:N);
for i = 1:N
    Sfi = Sf{i};
    for j = 1:k
        square_center_x = xp(i);
        square_center_y = Sfi(j);
        square_x = [square_center_x - 0.8*dx/2, square_center_x + 0.8*dx/2, square_center_x + 0.8*dx/2, square_center_x - 0.8*dx/2, square_center_x - 0.8*dx/2];
        square_y = [square_center_y - 1/2, square_center_y - 1/2, square_center_y + 1/2, square_center_y + 1/2, square_center_y - 1/2];
        if(i==1)&&(j==1)
            fill(square_x, square_y, 'g', 'FaceAlpha', '0.5', 'DisplayName','$S_f$'); % Fill the square
        else
            fill(square_x, square_y, 'g', 'FaceAlpha', '0.5', 'HandleVisibility','off'); % Fill the square
        end
    end
end

nexttile(tcl); hold all; legend('Location', 'southeast', 'interpreter','latex');
set(gca, 'FontName' , 'Times','FontSize',14); 
xlabel('$x$', 'Interpreter','latex', 'FontSize',24);
ylabel('$I_i$', 'Interpreter','latex', 'FontSize',24);
ylim([0.5,N+0.5]);
yticks(1:N);
for i = 1:N
    Sai = Sa{i};
    for j = 1:k
        square_center_x = xp(i);
        square_center_y = Sai(j);
        square_x = [square_center_x - 0.8*dx/2, square_center_x + 0.8*dx/2, square_center_x + 0.8*dx/2, square_center_x - 0.8*dx/2, square_center_x - 0.8*dx/2];
        square_y = [square_center_y - 1/2, square_center_y - 1/2, square_center_y + 1/2, square_center_y + 1/2, square_center_y - 1/2];
        if(i==1)&&(j==1)
            fill(square_x, square_y, 'b', 'FaceAlpha', '0.5', 'DisplayName','$S_a$'); % Fill the square
        else
            fill(square_x, square_y, 'b', 'FaceAlpha', '0.5', 'HandleVisibility','off'); % Fill the square
        end
    end
end