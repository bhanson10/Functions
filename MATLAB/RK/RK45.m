function [x,t,xf,tf]=RK45(f,x0,t0,s,p)
% RK45.M
% Benjamin Hanson, 2024
% 
% Simple implementation of the classical RK4/5 method
% 
% Inputs: 
%   f:   Equations of motion
%   x0:  Initial condition
%   t0:  Initial time
%   s:   Simulation parameters
%       - T:   Period
%       - dt:  Initial timestep
%       - eps: Error tolerance 
%   p:   Function parameters
%
% Outputs:
%   x:   Span of x
%   t:   Span of t
%   xf:  Final x
%   tf:  Final t

if ~isa(f, 'function_handle')
    error('Input must be a function handle (@f).');
end

if ~isfield(s,'dt')
    error('No initial timestep provided.')
end

if ~isfield(s,'T')
    error('No period provided.')
end

if ~isfield(s,'eps')
    error('No error tolerance provided.')
end

DT=s.dt; xf = x0; tf = t0; 
x = zeros(1,length(x0)); x(1,:) = xf; 
t = []; t(1) = tf;
while tf < s.T
    dt = DT/2; 
    f1=f(xf,p); 
    f2=f(xf+DT.*f1./2,p); 
    f3=f(xf+DT.*f2./2,p); 
    f4=f(xf+DT.*f3,p); 
    XF=xf+DT.*(f1./6+(f2+f3)./3+f4./6);  % calculate x using one RK4 step with timestep DT
    
    f1=f(xf,p); 
    f2=f(xf+dt.*f1./2,p); 
    f3=f(xf+dt.*f2./2,p); 
    f4=f(xf+dt.*f3,p); 
    xf=xf+dt.*(f1./6+(f2+f3)./3+f4./6);
    f1=f(xf,p); 
    f2=f(xf+dt.*f1./2,p); 
    f3=f(xf+dt.*f2./2,p); 
    f4=f(xf+dt.*f3,p); 
    xf=xf+dt.*(f1./6+(f2+f3)./3+f4./6);  % calculate x using two RK4 step with timestep DT/2
    
    delta = norm(xf-XF,1)/15;
    xf = (xf.*16-XF)./15; 
    tf = tf + DT;
    DT=min(DT*(DT*s.eps/delta)^(1/4),s.T-tf);

    x(end+1,:) = xf;
    t(end+1)   = tf;
end
end % function RK45