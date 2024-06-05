const.sigma=4; const.b=1; const.r=48; % Lorenz Parameters
x0  = [-11.5; -10; 9.5]; % [x;y;z]
s.dt = 0.01; s.T = 15; s.eps = 1e-5;

[x,t,xf,tf]=RK4(@Lorenz3D,x0,0,s,const); %[x,t,xf,tf]=RK45(@Lorenz3D,x0,0,s,const);

f2 = figure(2); clf; hold all; f2.Position = [750 200 600 475];
view(-109,14); lighting phong; light('Position',[-1 0 0]);
set(gca, 'FontName' , 'Times','FontSize',12);
xlabel("x", 'FontSize', 18, 'FontName', 'Times');
ylabel("y", 'FontSize', 18, 'FontName', 'Times');
zlabel("z", 'FontSize', 18, 'FontName', 'Times');
set(get(gca,'ZLabel'), 'Rotation', 0);
set(gca, 'FontName' , 'Times');
plot3(x(:,1),x(:,2),x(:,3),'Color', '[0 1 0 0.5]','linewidth',0.3, 'HandleVisibility','off');
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=Lorenz3D(y,const)                          
    f=[const.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -const.b*y(3)+y(1)*y(2)-const.b*const.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%