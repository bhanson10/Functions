function [xk,wk,xest,tspan] = PF(varargin)
% PF.M
% Benjamin Hanson, 2024
% 
% Perform Particle Filtering
% 
%   xest0:   Initial estimate (required)
%   Pest0:   Initial covariance (required)
%   dt:      Initial time step (required)
%   z:       Measurements wkth epochs (required)
%       - Measurements can be in the followkng two formats:
%           * z = {T} where T is the period, meaning no measurements
%           * z = {[t1,t2,t3;...,T],[z1,z2,z3,...,zn]} where ti are the
%           measurement epochs, and zi are the measurements
%   Q:       Process noise matrix (required)      
%   R:       Measurement noise matrix (optional if no measurements)      
%   f:       Dynamics model (required, function handle)
%   l:       Likelihood function (required, function handle)
%   h:       Measurement model (optional if no measurements, function handle)
%   method:  Time-marching method (optional)
%       - 'EE':    Explicit Euler
%       - 'RK4':   Runge-Kutta 4 (default)
%       - 'RK45':  Adaptive Runge-Kutta 4/5 
%   const:   Miscellaneous constants (optional)
%   num:     Number of particles (optional)
%   filter:  Filter type (SIS/SIR/Generic/...)
%       - If filter = 'SIS', Nthr = 0, resampling is not performed (default)
%       - If filter = 'SIR', Nthr = inf, resampling is performed every step 
%       - if filter = <number>, Nthr = number, resampling is performed
%       based on effective sample size 
%       - Else, Nthr is a threshold in [1,N]
%   particles: Pre-initialized particles (for different coordinate systems)
%   time:    Time frame of dynamics model, either "CT" or "DT" (required)
%
% Outputs:
%   xk:    Particles over time
%   wk:    Weights over time
%   tspan: Time span
%
% NOTE: RK45 does not work currently.

% Requirements
if ~(any(strcmp(varargin,'xest0'))&&any(strcmp(varargin,'Pest0'))&&any(strcmp(varargin,'dt'))&&any(strcmp(varargin,'z'))&&any(strcmp(varargin,'Q'))&&any(strcmp(varargin,'f'))&&any(strcmp(varargin,'l'))&&any(strcmp(varargin,'time')))
    error('Missing required components.')
end

% Default
R            = [];
h            = [];
pexist       = 0; 
const.exist  = 0;
method       = 'RK4';
num          = 500;
Nthr         = 0; 
particles    = 0; 

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
    elseif strcmp('l',varargin{i})
        l = varargin{i+1};
    elseif strcmp('h',varargin{i})
        h = varargin{i+1};
    elseif strcmp('method',varargin{i})
        method = varargin{i+1};
    elseif strcmp('const',varargin{i})
        const = varargin{i+1};
        const.exist = 1; 
    elseif strcmp('num',varargin{i})
        num = varargin{i+1};
    elseif strcmp('time',varargin{i})
        time = varargin{i+1};
    elseif strcmp('filter',varargin{i})
        filter = varargin{i+1};
    elseif strcmp('particles',varargin{i})
        particles = varargin{i+1};
        pexist = 1; 
    else
        error(append("Unspecified argument: ", varargin{i}));
    end
end

if strcmp(filter, 'SIS')
    Nthr = 0;
elseif strcmp(filter, 'SIR')
    Nthr = inf;
else
    Nthr = filter;
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
    error('Invalid time frame. Either "CT" or "DT".');
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

if ~isa(l, 'function_handle')
    error('Input must be a function handle (@l).');
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

% Initializing weights and particles
xk = {}; wk = {}; sum = 0; xest = [];
for k=1:num
    if ~pexist
        xk{end+1} = mvnrnd(xest0,Pest0)'; 
        wk{end+1} = l(xk{end}, xest0, Pest0, const);
    else
        xki = particles{1}; wki = particles{2};
        xk{end+1} = xki{k};
        wk{end+1} = wki{k};
    end
    sum = sum + wk{k};
end
xesti = zeros(size(xest0));
for k=1:num
    wj = wk{k};
    wk{k} = wj/sum;
    xesti = xesti + wk{k}.*xk{k};
end
xest(:,end+1) = xesti;

nz = size(z{2}); nz = nz(2);
t0 = 0; tspan = t0; 
for k=1:nz
    tk1 = z{1}; tk1 = tk1(k); zk1 = z{2}; zk1 = zk1(:,k); t0 = tspan(end);

    % Prediction
    dtx = dt;
    while t0 < tk1
        dtx = min(dtx, tk1-t0); t0 = t0 + dtx; tspan(end+1) = t0; xesti = zeros(size(xest0));
        for j=1:num
            xj = xk{j}; wj = wk{j};
            dfj = df(f,xj(:,end),dtx,const);
            xj(:,end+1) = mvnrnd(dfj{1},Q)'; 
            xk{j} = xj;
            wj(end+1) = wj(end);
            wk{j} = wj;
            xesti = xesti + xj(:,end).*wj(end);
        end
        xest(:,end+1) = xesti;
    end

    % Correction and Resampling
    if ~isnan(zk1)
        % Weight Correction
        sum = 0;
        for j=1:num
            xj = xk{j}; wj = wk{j};
            wj(end) = l(h(xj(:,end)), zk1, R, const);
            % if Nthr == inf
            %     wj(end) = l(h(xj(:,end)), zk1, R, const);
            % else
            %     wj(end) = wj(end)*l(h(xj(:,end)), zk1, R, const);
            % end
            wk{j} = wj; 
            sum = sum + wj(end);
        end
        xesti = zeros(size(xest0)); Neff = 0;
        for j=1:num
            xj = xk{j}; wj = wk{j};
            wj(end) = wj(end)/sum;
            Neff = Neff + wj(end)^2;
            wk{j} = wj; 
            xesti = xesti + xj(:,end).*wj(end);
        end
        Neff = 1/Neff; 
        xest(:,end) = xesti;
        
        % Resample based on ESS threshold
        if Neff < Nthr
            c = NaN(num,1);
            w1 = wk{1}; c(1) = w1(end);
            for j=2:num
                wj = wk{j};
                c(j) = c(j-1) + wj(end);
            end
            i = 1; 
            u1 = (1/num)*rand(1);
            ij = NaN(1,j);
            for j = 1:num
                uj = u1 + ((j-1)/num);
                ci = c(i);
                while uj > ci
                    i = i + 1;
                    ci = c(i);
                end
                xj = xk{j}; xki = xk{i}; xj(:,end) = xki(:,end); xk{j} = xj; 
                wj = wk{j}; wj(end) = (1/num); wk{j} = wj; 
                ij(j) = i; 
            end
        end
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