function oe = rv2oe(rv, mu, angletype, anomaly, oetype)
%    This function computes the classical (Keplerian) orbit
%    elements given the Cartesian position and velocity vectors, the
%    gravitational parameter of the attracting body, and the desired angletype
%
%    Parameters:
%        rv -- [Cartesian position vector [ km ] Cartesian velocity vector [ km/s ]]
%        mu -- gravitational parameter of attracting body [ km^3/s^2 ]
%        angletype --  string indicating desired angletype for inputs
%                      'deg' = degrees
%                      'rad' = radians
%        anomaly -- string indicating desired anomaly
%                      'f' = true anomaly
%                      'E' = eccentric anomaly 
%                      'M' = mean anomaly 
%        oetype -- string indiciating element type
%                      'coe' = classical (Keplerian) orbital elements (default)
%                      'eoe' = Equinoctal orbital elements
%
%    Outputs:
%        oe -- vector which contains the classical orbit elements
%               [ a; e; i; Om; w; anomaly ]
    
if ~exist('oetype', 'var'), oetype = 'coe'; end

r = rv(1:3);
v = rv(4:6);

oe = zeros(6,1); 
r2d = 180/pi();
I = [1 0 0];
J = [0 1 0];
K = [0 0 1]; 

rhat = r/norm(r); % position unit vector [km]
h = cross(r,v);   % angular momentum vector [km^2/s]
hhat = h/norm(h); % normalized angular momentum 
nhat = cross(K,h)/norm(cross(K,h)); % normalized ascending node vector

% Eccentricity
e = (1/mu)*cross(v,h)-rhat; % Eccentricity vector
oe(2) = norm(e); 

energy = (1/2)*dot(v,v)-(mu/norm(r)); % energy km^2/s^2
% If energy < 0, the orbit is closed (periodic)

% Semi-major axis (a) and parameter (p)
if(oe(2) ~= 1)
    oe(1) = -mu/(2*energy);
    p = oe(1)*(1-oe(2)^2);
else
    oe(1) = Inf;
    p = norm(h)^2/mu;
end

% Inclination (i) of orbit
oe(3) = acos(dot(K,hhat)); %If i < 90 deg, the elliptical orbit is a direct (prograde) orbit

% Right ascension of the ascending node (Omega)
oe(4) = mod(atan2(dot(J,nhat),dot(I,nhat)),2*pi()); 

% Argument of periapsis (w)
oe(5) = mod(atan2(dot(hhat,cross(nhat,e)),dot(nhat,e)),2*pi());

% True anomaly (f) at epoch [rad]
oe(6) = mod(atan2(norm(h)*dot(r,v),norm(h)^2-mu*norm(r)),2*pi());

if strcmp(oetype, 'coe')
    if~strcmp(anomaly,'f')
        % Eccentric anomaly (E) at epoch [rad]
        E = mod(2*atan2(sqrt(1-oe(2))*tan(oe(6)/2),sqrt(1+oe(2))),2*pi());
    
        if(strcmp(anomaly,'E'))
            oe(6) = E;
        elseif(strcmp(anomaly,'M'))
            % Mean anomaly (M) at epoch [rad]
            oe(6) = mod((E-oe(2)*sin(E)),2*pi()); 
        end
    end
    
    if(strcmp(angletype,'deg'))
        oe(3:6) = oe(3:6)*r2d;
    end
elseif strcmp(oetype, 'eoe')
    a = oe(1); e = oe(2); i = oe(3); Om = oe(4); om = oe(5); f = oe(6); 

    oe(1) = a * (1 - e^2);        % Semiparameter [ km ]
    oe(2) = e * cos(om + Om);     % x-component of eccentricity vector in orbital frame [ ] 
    oe(3) = e * sin(om + Om);     % y-component of eccentricity vector in orbital frame [ ]
    oe(4) = tan(i / 2) * cos(Om); % x-component of node vector in orbital frame [ ]
    oe(5) = tan(i / 2) * sin(Om); % y-component of node vector in orbital frame [ ]
    oe(6) = Om + om + f;          % True longitude [ rad ] 
    
    if(strcmp(angletype,'deg'))
        oe(6) = oe(6) * r2d;
    end
else
    error('Invalid orbital element type.')
end
end