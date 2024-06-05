clear all; close all; clc; 

%% Defining x_i and v_i 
v=@(x, const) (x<=const.mu)-(x>const.mu)/2; % Step Function
f=@(u) u.^2/2; % Burger's Equation

const.mu=0; 
xl=[-1, 1]; N=40; dx=(xl(2)-xl(1))/(N-1); xbar=xl(1):dx:xl(2);  
vbar=NaN(N,1); 
for i=1:N
    xim = xbar(i)-(dx/2); xip = xbar(i)+(dx/2);
    vbar(i) = (1/dx)*integral(@(x) v(x, const), xim, xip); % Cell Averages
end

%% Plotting
ttl.String = {'$t=$0'}; lbl.String = {'$x$', '$v$'}; lbl.FontSize = 24; 
lmts.YLimit = [-0.65,1.2];
initialize_figures('n', 1, 'margin', [800 200]', 'ttl', ttl, 'lgd', 1, 'lbl', lbl, 'lmts', lmts);

figure(1); 
p1 = plot(linspace(xl(1),xl(2),100000),v(linspace(xl(1),xl(2),100000),const),'r-','LineWidth',1,'DisplayName','$v(x)$');
p2 = scatter(xbar,vbar,200,'r','DisplayName','$\bar{v}_i$');
frames(1) = getframe(gcf); 

%% kth order ENO 1D FVM w/ Dirichlet boundary conditions
Tmax=0.5; dt=0.01; xh=xbar(1:N-1)+(dx/2); im=1:N-1; ip=2:N; 

k = 4; flux = 'Godunov'; t = 0; 
for i=1:Tmax/dt
    t=t+dt; const.mu=t; 
    V=ND1D(xbar,vbar,dx,k);
    [vm,Sm]=ENO1DNR(xbar,V,dx,xh,im,k); [vp,Sp]=ENO1DNR(xbar,V,dx,xh,ip,k); 
    if strcmp(flux,'Godunov')
        fh=Godunov_Flux_1D(f,vm,vp);
    end
    vbar(2:N-1)=vbar(2:N-1)-(dt/dx)*(fh(2:end)-fh(1:end-1));
    delete(p1); delete(p2);
    title(append('$t=$',num2str(t)),'Interpreter','latex', 'FontSize',24);
    p1 = plot(linspace(xl(1),xl(2),100000),v(linspace(xl(1),xl(2),100000),const),'r-','LineWidth',1,'DisplayName','$v(x)$');
    p2 = scatter(xbar,vbar,200,'r','DisplayName','$\bar{v}_i$');
    drawnow;
    frames(i+1) = getframe(gcf); 
end

create_video(frames, append('ENO1D_Burgers_',flux,'_',num2str(k),'.mp4'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = Godunov_Flux_1D(f,a,b)
    n=length(a); h=NaN(n,1); fn = @(x) -f(x); 
    for i=1:n
        if a(i)<=b(i)
            h(i)=fminbnd(f,a(i),b(i));
        else
            h(i)=fminbnd(fn,b(i),a(i));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_video(F, title)
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = 20;
    open(writerObj);
    for i=1:length(F)
        frame = F(i);
        writeVideo(writerObj, frame);    
    end
    close(writerObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%