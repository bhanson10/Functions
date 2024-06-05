function oe_cov =  oe2oe_cov(oe, oe_cov, anomaly, oetypein, oetypeout)
%    This function computes the uncertainty in the orbit elements 
%    given the Cartesian position and velocity vectors and uncertainty, the 
%    gravitational parameter of the attracting body, and the desired angletype
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
%        oe_cov    -- Orbital element covariance
%        angletype -- string indicating desired angletype for inputs
%                        'deg' = degrees, 
%                        'rad' = radians
%        anomaly   -- string indicating input anomaly
%                      'f' = true anomaly
%                      'E' = eccentric anomaly 
%                      'M' = mean anomaly 
%        oetypein  -- string indicating orbit element input type
%                        'coe' = classical (Keplerian)
%                        'eoe' = equinoctal
%        oetypeout -- string indicating orbit element input type
%                        'coe' = classical (Keplerian)
%                        'eoe' = equinoctal
%    Outputs:
%        oe_cov  -- Orbital element covariance

if strcmp(oetypein, oetypeout)
    return
else
    d2r = pi()/180; r2d = 180/pi();

    if strcmp(oetypein, 'coe')
        
        %  Classical orbit elements
        a_val  = oe(1); % semi-major axis [ km ]
        e_val  = oe(2); % eccentricity [ ]
        i_val  = oe(3); % inclination [ rad/deg ] 
        Om_val = oe(4); % longitude of the ascending node [ rad/deg ] 
        om_val = oe(5); % argument of periapsis [ rad/deg ] 
        A_val  = oe(6); % anomaly [ rad/deg ] 
        
        if(strcmp(anomaly,'f'))
            f_val = A_val;  
        elseif(strcmp(anomaly,'E'))
            f_val = 2*atan2(sqrt(1+e_val)*tan(A_val/2),sqrt(1-e_val)); 
        elseif(strcmp(anomaly,'M'))
            f_val = KepEqn(A_val, e_val, 'rad');
        end

        if strcmp(oetypeout, 'eoe')
            syms a e inc Om om nu 

            p = a*(1 - e^2);      % semiparameter [ km ] 
            f = e*cos(om + Om);   % x-component of eccentricity vector in orbital frame [ ]
            g = e*sin(om + Om);   % y-component of eccentricity vector in orbital frame [ ]
            h = tan(inc/2)*cos(Om); % x-component of node vector in orbital frame [ ]
            k = tan(inc/2)*sin(Om); % y-component of node vector in orbital frame [ ]
            L = Om + om + nu;     % true longitude [ rad ]

            J = jacobian([p,f,g,h,k,L],[a;e;inc;Om;om;nu]);

            oe_cov = double(subs(J,[a;e;inc;Om;om;nu],[a_val;e_val;i_val;Om_val;om_val;f_val])*oe_cov*subs(J',[a;e;inc;Om;om;nu],[a_val;e_val;i_val;Om_val;om_val;f_val])); 
        end
    elseif strcmp(oetypein, 'eoe')

        %  Equinoctal orbit elements
        p_val = oe(1); % semiparameter [ km ] 
        f_val = oe(2); % x-component of eccentricity vector in orbital frame [ ]
        g_val = oe(3); % y-component of eccentricity vector in orbital frame [ ]
        h_val = oe(4); % x-component of node vector in orbital frame [ ]
        k_val = oe(5); % y-component of node vector in orbital frame [ ]
        L_val = oe(6); % true longitude [ rad ]

        if strcmp(oetypeout, 'coe')
            syms p f g h k L
 
            a = p/(1 - f^2 - g^2);                       % semi-major axis [ km ]
            e = sqrt(f^2 + g^2);                         % eccentricity [ ]
            i = atan2(2*sqrt(h^2 + k^2), 1 - h^2 - k^2); % inclination [ rad ] 
            Om = atan2(k, h);                            % longitude of the ascending node [ rad ] 
            om = atan2(g*h - f*k, f*h + g*k);            % argument of periapsis [ rad ] 
            nu = L - atan2(g, f);                         % true anomaly [ rad ] 
            E = 2*atan2(sqrt(1-e)*tan(f/2),sqrt(1+e));   % eccentric anomaly [ rad ]
            M = E - e*sin(E);                            % mean anomaly [ rad ]
            
            if strcmp(anomaly, 'f')
                J = jacobian([a,e,i,Om,om,nu],[p;f;g;h;k;L]);
            elseif strcmp(anomaly, 'E')
                J = jacobian([a,e,i,Om,om,E],[p;f;g;h;k;L]);
            elseif strcmp(anomaly, 'M')
                J = jacobian([a,e,i,Om,om,M],[p;f;g;h;k;L]);
            end

            oe_cov = double(subs(J,[p;f;g;h;k;L],[p_val;f_val;g_val;h_val;k_val;L_val])*oe_cov*subs(J',[p;f;g;h;k;L],[p_val;f_val;g_val;h_val;k_val;L_val])); 
        end
    end
end
end