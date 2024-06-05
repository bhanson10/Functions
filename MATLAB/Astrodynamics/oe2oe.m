function oe = oe2oe(oe, angletype, anomaly, oetypein, oetypeout)
%    This function computes the orbit elements of one type
%    given the orbital elements of another type
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
%        angletype -- string indicating desired angletype for inputs
%                        'deg' = degrees, 
%                        'rad' = radians
%        anomaly -- string indicating input anomaly
%                      'f' = true anomaly
%                      'E' = eccentric anomaly 
%                      'M' = mean anomaly 
%        oetypein    -- string indicating orbit element input type
%                        'coe' = classical (Keplerian)
%                        'eoe' = equinoctal
%        oetypeout   -- string indicating orbit element input type
%                        'coe' = classical (Keplerian)
%                        'eoe' = equinoctal
%    Outputs:
%        rv     -- (Cartesian position vector [ km ], Cartesian velocity vector [ km/s ])  

if strcmp(oetypein, oetypeout)
    return
else
    d2r = pi()/180; r2d = 180/pi();
    
    if strcmp(oetypein,'coe')
        
        %  Classical orbit elements
        a  = oe(1); % semi-major axis [ km ]
        e  = oe(2); % eccentricity [ ]
        i  = oe(3); % inclination [ rad/deg ] 
        Om = oe(4); % longitude of the ascending node [ rad/deg ] 
        om = oe(5); % argument of periapsis [ rad/deg ] 
        A  = oe(6); % anomaly [ rad/deg ] 
        
        % Convert to radians if in degrees
        if(strcmp(angletype,'deg'))
            i   = i*d2r;
            Om  = Om*d2r;
            om  = om*d2r;
            A   = A*d2r;
        end
    
        if(strcmp(anomaly,'f'))
            f = A;  
        elseif(strcmp(anomaly,'E'))
            f = 2*atan2(sqrt(1+e)*tan(A/2),sqrt(1-e)); 
        elseif(strcmp(anomaly,'M'))
            f = KepEqn(A, e, 'rad');
        end
        
        %  Mod out by pi or 2pi
        i  = mod(i,pi());
        Om = mod(Om,(2*pi()));
        om = mod(om,(2*pi())); 
        f  = mod(f,(2*pi()));
        
        if strcmp(oetypeout,'eoe')
            oe(1) = a*(1 - e^2);      % semiparameter [ km ] 
            oe(2) = e*cos(om + Om);   % x-component of eccentricity vector in orbital frame [ ]
            oe(3) = e*sin(om + Om);   % y-component of eccentricity vector in orbital frame [ ]
            oe(4) = tan(i/2)*cos(Om); % x-component of node vector in orbital frame [ ]
            oe(5) = tan(i/2)*sin(Om); % y-component of node vector in orbital frame [ ]
            oe(6) = Om + om + f;      % true longitude [ rad ]
            
            % Convert to degrees
            if(strcmp(angletype,'deg'))
                oe(6) = oe(6)*r2d;
            end

        end
    elseif strcmp(oetypein,'eoe')
    
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
        L = mod(L,(2*pi()));


        if strcmp(oetypeout,'coe')
            oe(1) = p/(1 - f^2 - g^2);                       % semi-major axis [ km ]
            oe(2) = sqrt(f^2 + g^2);                         % eccentricity [ ]
            oe(3) = atan2(2*sqrt(h^2 + k^2), 1 - h^2 - k^2); % inclination [ rad ] 
            oe(4) = atan2(k, h);                             % longitude of the ascending node [ rad ] 
            oe(5) = atan2(g*h - f*k, f*h + g*k);             % argument of periapsis [ rad ] 
            oe(6) = L - atan2(g, f);                         % true anomaly [ rad ] 
    
            if(strcmp(anomaly,'E'))
                oe(6) = 2*atan2(sqrt(1-e)*tan(oe(6)/2),sqrt(1+e)); 
            elseif(strcmp(anomaly,'M'))
                oe(6) = 2*atan2(sqrt(1-e)*tan(oe(6)/2),sqrt(1+e)); 
                oe(6) = oe(6) - oe(2)*sin(oe(6));
            end
    
            % Convert to degrees
            if(strcmp(angletype,'deg'))
                oe(6) = oe(6)*r2d;
            end
        end
    end
end