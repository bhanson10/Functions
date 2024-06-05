function [x, z, tx, tz] = get_truth_measurements(varargin)
% get_measurements.m
% Benjamin Hanson, 2024
% 
% Given a dynamics model, initial state, and time parameters, return the
% true states and measurements over time.
% 
% Inputs: 
%   x0:      Initial true state (required)
%   dt:      True time step (required)
%   t:       Time parameter (required)
%       - Time parameter can be in the following two formats:
%           * t = T where T is period (no measurements)
%           * t = [dtz, T], where T is period and dtz is the measurement
%           time step
%           * t = [t1, t2, ..., T] where ti are the epochs
%   f:       Dynamics model (required, continuous-time, function handle)
%   h:       Measurement model (optional if no measurements, function handle)
%   Q:       Process noise matrix (required)
%   R:       Measurement noise matrix (optional if no measurements)
%   time:    Time frame of dynamics model, either "CT" or "DT" (required)
%   method:  Time-marching method (optional)
%       - 'EE':    Explicit Euler
%       - 'RK4':   Runge-Kutta 4 (default)
%       - 'RK45':  Adaptive Runge-Kutta 4/5
%   const:   Miscellaneous constants (optional)
%
% Outputs:
%   x:       True states
%   z:       Measurements
%   tx: True time span
%   tz: Measurement time span

% Requirements
if ~(any(strcmp(varargin,'x0'))&&any(strcmp(varargin,'dt'))&&any(strcmp(varargin,'t'))&&any(strcmp(varargin,'f'))&&any(strcmp(varargin,'Q'))&&any(strcmp(varargin,'time')))
    error('Missing required components.')
end

% Default
R           = [];
const.exist = 0;
method      = 'RK4';

% Optionals
for i=1:2:nargin
    if strcmp('x0',varargin{i})
        x0 = varargin{i+1};
    elseif strcmp('dt',varargin{i})
        dt = varargin{i+1};
    elseif strcmp('t',varargin{i})
        t = varargin{i+1};
    elseif strcmp('f',varargin{i})
        f = varargin{i+1};
    elseif strcmp('h',varargin{i})
        h = varargin{i+1};
    elseif strcmp('Q',varargin{i})
        Q = varargin{i+1};
    elseif strcmp('R',varargin{i})
        R = varargin{i+1};
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

if strcmp(time, 'CT')
    if strcmp(method, 'EE')
        df = @(f,x,dt,const) EE(f,x,dt,const);
    elseif strcmp(method, 'RK4')
        df = @(f,x,dt,const) RK4(f,x,dt,const);
    elseif strcmp(method, 'RK45')
        error('RK45 is not currently working.');
    else
        error('Invalid time-marching method.')
    end
elseif strcmp(time, 'DT')
    if ~const.exist
        df = @(f,x,dt,const) f(x,dt);
    else
        df = @(f,x,dt,const) f(x,dt,const);
    end
end

if length(t)==2
    dtz = t(1); T = t(2); 
    t = []; t0 = 0; 
    while (t0 < T)
        t0 = t0 + dtz; t(end+1) = t0; dtz = min(dtz,T-t0); 
    end
    if isempty(h)
        error('Missing measurement model.');
    end
elseif length(t)==1
    z = NaN;
else
    if isempty(h)
        error('Missing measurement model.');
    end
end

tx = 0; x = x0;
tz = t; nz = length(tz); z = [];

t0 = 0; 
for i=1:nz
    tk1 = tz(i); 
    
    dtx = dt;
    % True State
    while t0 < tk1
        dtx = min(dtx, tk1-t0); t0 = t0 + dtx; tx(end+1) = t0;
        w = sqrt(Q)*randn(length(Q),1); 
        x(:,end+1) = df(f,x(:,end),dtx,const) + w;
    end
    
    if length(t)~=1
        % Measurement
        v = sqrt(R)*randn(length(R),1);
        z(:,end+1) = h(x(:,end)) + v;
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
function x1 = RK4(f,x,dt,const)
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%