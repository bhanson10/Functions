function [x,t,xf,tf]=RK4(f,x0,t0,s,p)
% RK4.M
% Benjamin Hanson, 2024
% 
% Simple implementation of the classical RK4 method
% 
% Inputs: 
%   f:   Equations of motion
%   x0:  Initial condition
%   t0:  Initial time
%   s:   Simulation parameters
%       - T:   Period
%       - dt:  Timestep
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
    error('No timestep selected.')
end

if ~isfield(s,'T')
    error('No period selected.')
end

dt=s.dt; xf = x0; tf = t0; 
x = zeros(1,length(x0)); x(1,:) = xf; 
t = []; t(1) = tf;
while tf < s.T
    f1=f(xf,p); 
    f2=f(xf+dt.*f1./2,p); 
    f3=f(xf+dt.*f2./2,p); 
    f4=f(xf+dt.*f3,p); 
    xf=xf+dt.*(f1./6+(f2+f3)./3+f4./6); 
    tf=tf+dt;
    dt = min(dt, s.T-tf);

    x(end+1,:) = xf;
    t(end+1)   = tf;
end
end % function RK4