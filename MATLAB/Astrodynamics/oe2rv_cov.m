function rv_cov =  oe2rv_cov(oe, oe_cov, mu_val, anomaly, oetype)
%    This function computes the uncertainty in the orbit elements 
%    given the Cartesian position and velocity vectors and uncertainty, the 
%    gravitational parameter of the attracting body, and the desired angletype
%
%    Parameters:
%        oe      -- Orbital elements [a [ km ]; e [ ]; i [ rad ]; Om [ rad ]; om [ rad ]; A [ rad ]];
%        oe_cov  -- Orbital element covariance
%        mu      -- gravitational parameter of attracting body [ km^3/s^2 ]
%        anomaly -- string indicating input anomaly
%                      'f' = true anomaly
%                      'E' = eccentric anomaly 
%                      'M' = mean anomaly 
%        oetype  -- string indicating orbit element input type
%                      'coe' = classical (Keplerian)
%                      'eoe' = equinoctal
%    Outputs:
%        rv_cov  -- Orbital element covariance

if ~exist('oetype', 'var'), oetype = 'coe'; end

if strcmp(oetype, 'eoe')
    syms p f g h k L

    a  = p / (1 - f^2  - g^2);
    e  = sqrt(f^2 + g^2);
    i  = atan2(2 * sqrt(h^2 + k^2), 1 - h^2 - k^2); 
    Om = atan2(k, h); 
    om = atan2(g * h - f * k, f * h + g * k); 
    nu  = L - atan2(g, f); 

    J = jacobian([a, e, i, Om, om, nu],[p; f; g; h; k; L]);
    
    oe_cov = double(subs(J,[p; f; g; h; k; L],[oe(1); oe(2); oe(3); oe(4); oe(5); oe(6)])*oe_cov*subs(J',[p; f; g; h; k; L],[oe(1); oe(2); oe(3); oe(4); oe(5); oe(6)]));

    a_val  = oe(1) / (1 - oe(2)^2  - oe(3)^2);
    e_val  = sqrt(oe(2)^2 + oe(3)^2);
    inc_val  = atan2(2 * sqrt(oe(4)^2 + oe(5)^2), 1 - oe(4)^2 - oe(5)^2); 
    Om_val = atan2(oe(5), oe(4)); 
    om_val = atan2(oe(3) * oe(4) - oe(2) * oe(5), oe(2) * oe(4) + oe(3) * oe(5)); 
    A_val = oe(6) - atan2(oe(3), oe(2)); anomaly = 'f'; 

    clear p f g h k L a e i Om om nu J 
elseif strcmp(oetype, 'coe')
    a_val = oe(1); e_val = oe(2); inc_val = oe(3); Om_val = oe(4); om_val = oe(5); A_val = oe(6);
end

syms mu a e inc Om om f E real

ehat  = [cos(om)*cos(Om)-cos(inc)*sin(om)*sin(Om); 
        cos(om)*sin(Om)+cos(inc)*sin(om)*cos(Om);
        sin(om)*sin(inc)];
ehatp = [-(sin(om)*cos(Om)+cos(inc)*cos(om)*sin(Om));
         -(sin(om)*sin(Om)-cos(inc)*cos(om)*cos(Om));
         cos(om)*sin(inc)];

if(strcmp(anomaly,'f'))
    r = (a*(1-e^2)/(1+e*cos(f))).*(cos(f).*ehat + sin(f).*ehatp); 
    v = sqrt(mu/(a*(1-e^2))).*(-sin(f).*ehat + (e+cos(f)).*ehatp);
    J = jacobian([r;v],[a;e;inc;Om;om;f]); 
    rv_cov = double(subs(J,[a;e;inc;Om;om;f;mu],[a_val;e_val;inc_val;Om_val;om_val;A_val;mu_val])*oe_cov*subs(J',[a;e;inc;Om;om;f;mu],[a_val;e_val;inc_val;Om_val;om_val;A_val;mu_val]));
elseif(strcmp(anomaly,'E'))
    f = 2*atan2(sqrt(1+e)*tan(E/2), sqrt(1-e));
    r = (a*(1-e^2)/(1+e*cos(f))).*(cos(f).*ehat + sin(f).*ehatp); 
    v = sqrt(mu/(a*(1-e^2))).*(-sin(f).*ehat + (e+cos(f)).*ehatp);
    J = jacobian([r;v],[a;e;inc;Om;om;E]); 
    rv_cov = double(subs(J,[a;e;inc;Om;om;E;mu],[a_val;e_val;inc_val;Om_val;om_val;A_val;mu_val])*oe_cov*subs(J',[a;e;inc;Om;om;E;mu],[a_val;e_val;inc_val;Om_val;om_val;A_val;mu_val]));
elseif(strcmp(anomaly,'M'))
    % Because there is not a close form for E in terms of M, this
    % Jacobian must be approximated
    
    n = numel(oe);
    m = numel(oe2rv(oe,mu_val,'rad','M','coe'));
    eps = 1e-8; 

    J = zeros(m,n);
    for count = 1:n
        oe_perturbed = oe; 
        oe_perturbed(count) = oe_perturbed(count) + eps; 
        delta_f = oe2rv(oe_perturbed,mu_val,'rad','M','coe') - oe2rv(oe,mu_val,'rad','M','coe');
        J(:,count) = delta_f./eps; 
    end
    rv_cov = J*oe_cov*J';
end
end
