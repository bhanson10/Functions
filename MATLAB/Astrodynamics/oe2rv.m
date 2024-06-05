function rv = oe2rv(oe, mu, angletype, anomaly, oetype)
%    This function computes the orbit elements 
%    given the Cartesian position and velocity vectors, the gravitational
%    parameter of the attracting body, and the desired angletype
%
%    Parameters:
%        oe        -- vector which contains the classical orbit elements
%                        [ a; e; i; Om; om; A ]
%                            a  - semimajor axis [ km ]
%                            e  - eccentricity [ ]
%                            i  - inclination [ rad/deg ]
%                            Om - longitude of the ascending node [ rad/deg ]
%                            om - argument of periapsis [ rad/deg ]
%                            A  - anomaly at epoch [ rad/deg ]
%                    OR the equinoctal orbit elements
%                        [ p; f; g; h; k; L ]
%                            p - semiparameter [ km ]
%                            f - x-component of the eccentricity vector in orbital frame [ ]
%                            g - y-component of the eccentricity vector in orbital frame [ ]
%                            h - x-component of the node vector in orbital frame [ ]
%                            k - y-component of the node vector in orbital frame [ ]
%                            L - true longitude [ rad/deg ]
%        mu        -- gravitational parameter of attracting body [ km^3/s^2 ]
%        angletype -- string indicating desired angletype for inputs
%                        'deg' = degrees, 
%                        'rad' = radians
%        anomaly -- string indicating input anomaly
%                      'f' = true anomaly
%                      'E' = eccentric anomaly 
%                      'M' = mean anomaly 
%        oetype    -- string indicating orbit element input type
%                        'coe' = classical (Keplerian)
%                        'eoe' = equinoctal
%    Outputs:
%        rv     -- (Cartesian position vector [ km ], Cartesian velocity vector [ km/s ])  

d2r = pi()/180; 

if strcmp(oetype,'coe')
    
    %  Classical orbit elements
    a  = oe(1); % semi-major axis
    e  = oe(2); % eccentricity
    i  = oe(3); % inclination
    Om = oe(4); % longitude of the ascending node
    om = oe(5); % argument of periapsis
    A  = oe(6); % anomaly
    
    % Convert to radians if in degrees
    if(strcmp(angletype,'deg'))
        i   = i*d2r;
        Om  = Om*d2r;
        om  = om*d2r;
        A   = A*d2r;
    end

    if(strcmp(anomaly,'f'))
        E = 2*atan2(sqrt(1-e)*tan(A/2),sqrt(1+e));
        M = E - e*sin(E); 
    elseif(strcmp(anomaly,'E'))
        M = A - e*sin(A); 
    elseif(strcmp(anomaly,'M'))
        M = A; 
    end
    
    %  Mod out by pi or 2pi
    i  = mod(i,pi());
    Om = mod(Om,(2*pi()));
    om = mod(om,(2*pi())); 
    M  = mod(M,(2*pi()));
    
    if (a > 0)  % ----- eccentric orbit ----- %
        % True anaomaly [ rad ]
        f = KepEqn(M, e, 'rad');
    else        % ----- hyperbolic orbit ----- %
        
        %  Compute hyperbolic anomaly
        %  A priori estimate
        j = 1;
        H = [M];
    
        %  Newton iteration to find hyperbolic anomaly
        %  Algorithm [goal: find H so f = 0]
        f_H = [e*sinh(H(j)) - H(j) - M];
        while abs(f_H(j)) > 1e-13
            H(end+1) = (H(j) - f_H(j)/(e*cosh(H(j)) - 1));
            j = j + 1;
            f_H(end+1) = (e*sinh(H(j)) - H(j) - M);
        end
    
        % Converged eccentric anomaly [ rad ]
        H = H(j);
    
        %  True anomaly [ rad ]
        f = mod((2*atan(sqrt((e + 1)/(e - 1))*tanh(H/2))),(2*pi()));
    end
    
    % Argument of latitude [ rad ]
    th = om + f; 
    
    % Magnitude of position vector [ km ]
    r = a*(1-e^2)/(1+e*cos(f)); % trajectory equation
    
    % Magnitude of velocity vector [ km/s ]
    v = sqrt(mu*(2/r-1/a)); % vis-viva equation
    
    % Ascending node vector
    nhat = [cos(Om); sin(Om); 0];
    rT   = [-cos(i)*sin(Om); cos(i)*cos(Om); sin(i)]; 
    
    gamma = atan2(e*sin(f), 1+e*cos(f)); % [ rad ]
    
    % Normalized position and velocity vectors
    rhat = cos(th)*nhat + sin(th)*rT; % [ km ]
    vhat = sin(gamma - th)*nhat + cos(gamma - th)*rT; % [ km/s ]
    
    % Position and velocity vectors
    r = r*rhat; % [ km ]
    v = v*vhat; % [ km/s ]
    
    rv = [r; v]; 

elseif strcmp(oetype,'eoe')

    % Equinoctal orbit elements
    p = oe(1); % semiparameter [ km ]
    f = oe(2); % x-component of eccentricity vector in orbital frame [ ]
    g = oe(3); % y-component of eccentricity vector in orbital frame [ ]
    h = oe(4); % x-component of node vector in orbital frame [ ]
    k = oe(5); % y-component of node vector in orbital frame [ ]
    L = oe(6); % true longitude [ rad/deg ]

    % Convert to radians if in degrees
    if(strcmp(angletype,'deg'))
        L   = L*d2r;
    end
    
    %  Mod out by 2pi 
    L  = mod(L,(2*pi()));

    alpha_sq = h^2 - k^2; 
    s_sq     = 1 + h^2 + k^2; 
    w        = 1 + f*cos(L) + g*sin(L); 
    r        = p/w; 

    r = [(r/s_sq)*(cos(L) + alpha_sq*cos(L) + 2*h*k*sin(L));...
         (r/s_sq)*(sin(L) - alpha_sq*sin(L) + 2*h*k*cos(L));...
         (2*r/s_sq)*(h*sin(L) - k*cos(L))];
    v = [(-1/s_sq)*sqrt(mu/p)*(sin(L) + alpha_sq*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha_sq*g);...
         (-1/s_sq)*sqrt(mu/p)*(-cos(L) + alpha_sq*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha_sq*f);...
         (2/s_sq)*sqrt(mu/p)*(h*cos(L) + k*sin(L) + f*h + g*k)]; 

    rv = [r; v]; 
end
end