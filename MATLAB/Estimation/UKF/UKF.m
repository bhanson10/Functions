function [xest,Pest,tspan] = UKF(varargin)
% UKF.M
% Benjamin Hanson, 2024
% 
% Perform Unscented Kalman Filter estimation given by Wan et. al. (2000)
% 
%   xest0:   Initial estimate (required)
%   Pest0:   Initial covariance (required)
%   dt:      Initial time step (required)
%   z:       Measurements with epochs (required)
%       - Measurements can be in the following two formats:
%           * z = {T} where T is the period, meaning no measurements
%           * z = {[t1,t2,t3;...,T],[z1,z2,z3,...,zn]} where ti are the
%           measurement epochs, and zi are the measurements
%   Q:       Process noise matrix (required)
%   R:       Measurement noise matrix (optional if no measurements)      
%   f:       Dynamics model (required, function handle)
%   time:    Time frame of dynamics model, either "CT" or "DT" (required)
%   h:       Measurement model (optional if no measurements, function handle)
%   method:  Time-marching method (optional)
%       - 'EE':    Explicit Euler 
%       - 'RK4':   Runge-Kutta 4 (default)
%       - 'RK45':  Adaptive Runge-Kutta 4/5 
%   const:   Miscellaneous constants (optional)
%   alpha:   Alpha, spread of the sigma points around xest0 (optional)
%   beta:    Beta, incorporates prior knowledge (optional)
%   kappa:   Kappa, secondary scaling parameter (optional)
%
% Outputs:
%   xest:    True states
%   Pest:    Estimated covariances
%   tspan_x: Time span
%
% NOTE: RK45 does not work currently.

% Requirements
if ~(any(strcmp(varargin,'xest0'))&&any(strcmp(varargin,'Pest0'))&&any(strcmp(varargin,'dt'))&&any(strcmp(varargin,'z'))&&any(strcmp(varargin,'Q'))&&any(strcmp(varargin,'f'))&&any(strcmp(varargin,'time')))
    error('Missing required components.')
end

% Default
R           = [];
h           = [];
const.exist = 0;
method      = 'RK4';
alpha       = 1e-3;
beta        = 2;
kappa       = 0;

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
    elseif strcmp('h',varargin{i})
        h = varargin{i+1};
    elseif strcmp('method',varargin{i})
        method = varargin{i+1};
    elseif strcmp('const',varargin{i})
        const = varargin{i+1};
        const.exist = 1; 
    elseif strcmp('alpha',varargin{i})
        alpha = varargin{i+1};
    elseif strcmp('beta',varargin{i})
        beta = varargin{i+1};
    elseif strcmp('kappa',varargin{i})
        kappa = varargin{i+1};
    elseif strcmp('time',varargin{i})
        time = varargin{i+1};
    else
        error(append("Unspecified argument: ", varargin{i}));
    end
end

if strcmp(time, 'CT')
    if strcmp(method, 'EE')
        df = @(f,x,dt,const) EE(f,x,dt,const);
    elseif strcmp(method, 'RK4')
        df = @(f,x,dt,const) RK4(f,x,dt,const);
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
else
    error('Invalid time frame. Either "CT" or "DT".')
end

% Regulations 
[rows, columns] = size(Q);
if ~(rows==columns) 
    error('Process noise matrix is not square.'); 
end

[rows, columns] = size(R);
if ~(rows==columns) 
    error('Measurement noise matrix is not square.'); 
end

[rows, columns] = size(Pest0);
if ~(rows==columns) 
    error('Covariance matrix is not square.'); 
end

if ~(issymmetric(Pest0))
    error('Covariance matrix is not symmetric.')
end

if (kappa+length(xest0)==0)
    error('Scaling parameter (k) plus dimension of xest0 (n_a) cannot equal 0.');
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

% Unscented Transform - Weights and Sigma Points
n   = length(xest0);          % Dimensions of xest0
lam = alpha^2*(n+kappa)-n; % Lambda
Wm  = NaN(1,2*n+1);        % Weights, m
Wc  = NaN(1,2*n+1);        % Weights, c

% Defining Weights/Sigma Points
Wm(1) = lam/(n + lam); 
Wc(1) = (lam/(n + lam)) + (1-alpha^2+beta); 
for i=2:2*n+1
    Wm(i) = 1/(2*(n+lam)); 
    Wc(i) = Wm(i);  
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

        % Calculate Sigma Points
        A = NaN(n,2*n+1);
        L = sqrtm((n+lam).*Pest{end});
        A(:,1) = xest(:,end); 
        for j=2:(n+1)
            A(:,j) = xest(:,end) + L(j-1,:)'; 
        end
        for j=(n+2):(2*n+1)
            A(:,j) = xest(:,end) - L(j-n-1,:)'; 
        end
        
        % Prediction
        for j = 1:2*n+1
            dfi = df(f,A(:,j),dt,const);
            A(:,j) = dfi{1}; % Sigma Points
        end
    
        xest(:,end+1) = zeros(length(xest0),1); 
        for j = 1:2*n+1
            xest(:,end) = xest(:,end) + Wm(j).*A(:,j); % Mean
        end
        
        Pest{end+1} = Q;
        for j = 1:2*n+1
            Pest{end} = Pest{end} + Wc(j).*((A(:,j)-xest(:,end))*(A(:,j)-xest(:,end))'); % Covariance
        end
    end

    % Correction
    if ~isnan(zk1)
        for j = 1:2*n+1
            Z(:,j) = h(A(:,j)); % Measurement
        end
    
        zpred = zeros(length(zk1),1);
        for j = 1:2*n+1
            zpred = zpred + Wm(j).*Z(:,j); % Measurement
        end
    
        P_zz = zeros(length(zk1));
        for j = 1:2*n+1
            P_zz = P_zz + Wc(j).*((Z(:,j)-zpred)*(Z(:,j)-zpred)');
        end
    
        P_xz = zeros(length(xest0),length(zk1));
        for j = 1:2*n+1
            P_xz = P_xz + Wc(j).*((A(:,j)-xest(:,end))*(Z(:,j)-zpred)');
        end
    
        S = R + P_zz; 
        K = P_xz*(S^(-1));
        xest(:,end) = xest(:,end) + K*(zk1-zpred);
        Pest{end}   = Pest{end} - K*S*K';

        
    end
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