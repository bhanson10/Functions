function oe_cov =  rv2oe_cov(rv, rv_cov, mu_val, anomaly, oetype, order)
%    This function computes the uncertainty in the orbit elements 
%    given the Cartesian position and velocity vectors and uncertainty, the 
%    gravitational parameter of the attracting body, and the desired angletype
%
%    Parameters:
%        rv      -- [Cartesian position vector [ km ]; Cartesian velocity vector [ km/s ]]
%        rv_cov  -- Cartesian covariance
%        mu      -- gravitational parameter of attracting body [ km^3/s^2 ]
%        anomaly -- string indicating input anomaly
%                      'f' = true anomaly
%                      'E' = eccentric anomaly 
%                      'M' = mean anomaly 
%        oetype  -- string indicating orbit element input type
%                      'coe' = classical (Keplerian)
%                      'eoe' = equinoctal
%        order   -- order of approximation, 1 or 2
%
%    Outputs:
%        oe_cov  -- Orbital element covariance

if ~exist('oetype', 'var'), oetype = 'coe'; end
if ~exist('order', 'var'), order = 1; end

r_val = rv(1:3); v_val = rv(4:6); 

syms mu real; r = sym('r',[3,1],'real'); v = sym('v',[3,1],'real');

xhat = [1; 0; 0]; yhat = [0; 1; 0]; zhat = [0; 0; 1]; % unit vectors

rhat = r./norm(r);                           % Unit vector r
En   = (1/2)*dot(v,v)-(mu/(dot(r,r)^(1/2))); % Energy
H    = cross(r,v);                           % Angular Momentum
hhat = H./norm(H);                           % Unit vector h

a    = -(mu/(2*En));                                            % Semi-major axis
evec = (1/mu).*cross(v,H)-rhat; e = norm(evec); ehat = evec./e; % Eccentricity
i  = acos(dot(hhat,zhat));                                      % Inclination
Om   = atan2(dot(hhat,xhat),-dot(hhat,yhat));                   % Long. of the asc. node
om   = atan2(dot(ehat,zhat),dot(cross(hhat,ehat),zhat));        % Arg. of periapsis
nu    = atan2(dot(hhat,cross(ehat,r)),dot(ehat,r));             % True anomaly
E    = 2*atan2(sqrt(1-e)*tan(nu/2),sqrt(1+e));                  % Eccentric anomaly
M    = E - e*sin(E);                                            % Mean anomaly

if strcmp(oetype, 'coe')
    if strcmp(anomaly, 'f')
        J = jacobian([a,e,i,Om,om,nu],[r;v]);
    elseif strcmp(anomaly, 'E')
        J = jacobian([a,e,i,Om,om,E],[r;v]);
    elseif strcmp(anomaly, 'M')
        J = jacobian([a,e,i,Om,om,M],[r;v]);
    end
elseif strcmp(oetype, 'eoe')
    p = a * (1 - e^2);        % Semiparameter
    f = e * cos(om + Om);     % x-component of eccentricity vector in orbital frame
    g = e * sin(om + Om);     % y-component of eccentricity vector in orbital frame
    h = tan(i / 2) * cos(Om); % x-component of node vector in orbital frame
    k = tan(i / 2) * sin(Om); % y-component of node vector in orbital frame
    L = Om + om + nu;          % True longitude 

    J = jacobian([p,f,g,h,k,L],[r;v]);
end

if order==1
    oe_cov = double(subs(J,[r;v;mu],[r_val;v_val;mu_val])*rv_cov*subs(J',[r;v;mu],[r_val;v_val;mu_val])); 
end
end