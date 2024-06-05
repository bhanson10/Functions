function [xest,Pest,tspan] = EKF(varargin)
% EKF.M
% Benjamin Hanson, 2024
% 
% Perform Extended Kalman Filter estimation given dynamics and measurement
% models
% 
% Inputs: 
%   xest0:   Initial estimate (required)
%   Pest0:   Initial covariance (required)
%   dt:      Initial time step (required)
%   z:       Measurements with epochs (required)
%       - Measurements can be in the following two formats:
%           * z = {T} where T is the period, meaning no measurements
%           * z = {[t1,t2,t3,...,T],[z1,z2,z3,...,zn]} where ti are the
%           measurement epochs, and zi are the measurements
%   Q:       Process noise matrix (required)
%   R:       Measurement noise matrix (optional if no measurements)
%   f:       Dynamics model (required, function handle)
%   time:    Time frame of dynamics model, either "CT" or "DT" (required)
%   h:       Measurement model (optional if no measurements, function handle)
%   F:       Dynamics Jacobian (optional, function handle)
%   H:       Measurement Jacobian (optional, function handle)
%   method:  Time-marching method (optional)
%       - 'EE':    Explicit Euler - Dynamics Jacobian is used
%       - 'RK4':   Runge-Kutta 4 (default) - Dynamics Jacobian is estimated
%       - 'RK45':  Adaptive Runge-Kutta 4/5 - Dynamics Jacobian is
%       estimated
%   const:   Miscellaneous constants (optional)

% Outputs:
%   xest:    Estimated states
%   Pest:    Estimated covariances
%   tspan:   Time span
%
% NOTE: RK45 does not work currently.

% Requirements
if ~(any(strcmp(varargin,'xest0'))&&any(strcmp(varargin,'Pest0'))&&any(strcmp(varargin,'dt'))&&any(strcmp(varargin,'z'))&&any(strcmp(varargin,'Q'))&&any(strcmp(varargin,'f'))&&any(strcmp(varargin,'time')))
    error('Missing required components.')
end

% Default
F           = [];
R           = [];
h           = [];
H           = [];
const.exist = 0; 
method      = 'RK4';

% Optionals
for i=1:2:nargin
    if strcmp('xest0',varargin{i})
        xest0 = varargin{i+1};
    elseif strcmp('Pest0',varargin{i})
        Pest0 = varargin{i+1};
    elseif strcmp('dt',varargin{i})
        dt = varargin{i+1};
    elseif strcmp('z',varargin{i})
        z = varargin{i+1};
    elseif strcmp('Q',varargin{i})
        Q = varargin{i+1};
    elseif strcmp('R',varargin{i})
        R = varargin{i+1};
    elseif strcmp('f',varargin{i})
        f = varargin{i+1};
    elseif strcmp('F',varargin{i})
        F = varargin{i+1};
    elseif strcmp('h',varargin{i})
        h = varargin{i+1};
    elseif strcmp('H',varargin{i})
        H = varargin{i+1};
    elseif strcmp('method',varargin{i})
        method = varargin{i+1};
    elseif strcmp('const',varargin{i})
        const = varargin{i+1};
        const.exist = 1; 
    elseif strcmp('time',varargin{i})
        time = varargin{i+1};
    else
        error(append("Unspecified argument: ", varargin{i}));
    end
end

% Checks and Balances
[rows, columns] = size(Q);
if ~(rows==columns) 
    error('Process noise matrix is not square.'); 
end

if ~(length(xest0)==rows)
    error('State vector and process noise matrix are not compatible dimensions.');
end

[rows, columns] = size(R);
if ~(rows==columns) 
    error('Measurement noise matrix is not square.'); 
end

[rows, columns] = size(Pest0);
if ~(rows==columns) 
    error('Covariance matrix is not square.'); 
end

if ~(length(xest0)==rows)
    error('State vector and covariance matrix are not compatible dimensions.');
end

if ~(issymmetric(Pest0))
    error('Covariance matrix is not symmetric.')
end

if strcmp(time, 'CT')
    if strcmp(method, 'EE')
        df = @(f,x,dt,const) EE(f,x,dt,const);
        if ~isempty(F)
            if ~const.exist
                dF = @(f,df,F,x,dt,const) eye(length(x))+dt.*F(x);
            else
                dF = @(f,df,F,x,dt,const) eye(length(x))+dt.*F(x,const);
            end
        else
            dF = @(f,df,F,x,dt,const) est_jacob_F(f,df,F,x,dt,const);  
        end
    elseif strcmp(method, 'RK4')
        df = @(f,x,dt,const) RK4(f,x,dt,const);
        dF = @(f,df,F,x,dt,const) est_jacob_F(f,df,F,x,dt,const); 
    elseif strcmp(method, 'RK45')
        error('RK45 not currently working.')
    else
        error('Invalid time-marching method.')
    end
elseif strcmp(time, 'DT')
    if ~const.exist
        df = @(f,x,dt,const) f(x,dt);
    else
        df = @(f,x,dt,const) f(x,dt,const);
    end
    if ~isempty(F)
        if ~const.exist
            dF = @(f,df,F,x,dt,const) F(x,dt);
        else
            dF = @(f,df,F,x,dt,const) F(x,dt,const);
        end
    else
        dF = @(f,df,F,x,dt,const) est_jacob_F(f,df,F,x,dt,const);  
    end
else
    error('Invalid time frame. Either "CT" or "DT".')
end

if isempty(H)&&~isempty(h)
    H = @(h,x) est_jacob_H(h,x);
else
    H = @(h,x) H(x);
end

if length(z)==1
    z = {z{1},NaN};
else
    if isempty(h)
        error("Included measurements but missing measurement model.")
    end
    if isempty(R)
        error("Included measurements but missing measurement noise.")
    end
end

xest = xest0;
Pest = {}; Pest{1} = Pest0;
nz = size(z{2}); nz = nz(2);
t0 = 0; tspan = t0;

for i=1:nz
    tk1 = z{1}; tk1 = tk1(i); zk1 = z{2}; zk1 = zk1(:,i); t0 = tspan(end);
    
    % Prediction
    dtx = dt;
    while t0 < tk1
        dtx = min(dtx, tk1-t0); t0 = t0 + dtx; tspan(end+1) = t0;
        dfi = df(f,xest(:,end),dtx,const);
        xest(:,end+1) = dfi{1}; dtx = dfi{2};
        Pest{end+1}   = dF(f,df,F,xest(:,end),dtx,const)*Pest{end}*dF(f,df,F,xest(:,end),dtx,const)'+Q;
    end

    % Correction
    if ~isnan(zk1)
        zpred = h(xest(:,end));
        K = Pest{end}*H(h,xest(:,end))'*(H(h,xest(:,end))*Pest{end}*H(h,xest(:,end))' + R)^(-1);
        xest(:,end) = xest(:,end) + K*(zk1-zpred);
        Pest{end} = (eye(length(xest(:,end)))-K*H(h,xest(:,end)))*Pest{end}; 
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dF = est_jacob_F(f,df,F,x,dt,const)
    dx = 1e-8;
    n = length(x);
    df0 = df(f,x,dt,const); df0 = df0{1};
    dF = zeros(n);
    for i=1:n
        if(x(i)~=0), Dxj = abs(x(i))*dx; else, Dxj = dx; end
        x_plus = zeros(size(x));
        for j=1:n
            if i~=j, x_plus(j) = x(j); else, x_plus(j) = x(j)+Dxj; end
        end
        dfi = df(f,x_plus,dt,const); dfi = dfi{1};
        dF(:,i) = (dfi-df0)./Dxj;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = est_jacob_H(h,x)
    dx = 1e-8;
    n = length(x);
    h0 = h(x);
    H = zeros(length(h0),n);
    for i=1:n
        if(x(i)~=0), Dxj = abs(x(i))*dx; else, Dxj = dx; end
        x_plus = zeros(size(x));
        for j=1:n
            if i~=j, x_plus(j) = x(j); else, x_plus(j) = x(j)+Dxj; end
        end
        H(:,i) = (h(x_plus)-h0)./Dxj;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = EE(f,x,dt,const)
    if ~isa(f, 'function_handle')
        error('Input must be a function handle (@f).');
    end
    
    if ~const.exist
        x1 = x + dt*f(x);
    else
        x1 = x + dt*f(x,const);
    end
    df = {x1,dt};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = RK4(f,x,dt,const)
    if ~isa(f, 'function_handle')
        error('Input must be a function handle (@f).');
    end
    
    if ~const.exist
        f1 = f(x);
        f2 = f(x+(dt/2).*f1);
        f3 = f(x+(dt/2).*f2);
        f4 = f(x+dt.*f3);
    else 
        f1 = f(x,const);
        f2 = f(x+(dt/2).*f1,const);
        f3 = f(x+(dt/2).*f2,const);
        f4 = f(x+dt.*f3,const);
    end
    
    x1 = x + dt.*((f1./6)+(f2./3)+(f3./3)+(f4./6));
    df = {x1,dt};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%