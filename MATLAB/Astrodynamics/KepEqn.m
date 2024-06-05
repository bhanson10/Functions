function f = KepEqn(M, e, angletype)
%Solves Kepler's equation 
%
%    Parameters:
%    M          -- mean anomaly at epoch 
%    e          -- eccentricity
%    angletype  -- string indicating desired angletype for inputs
%                    'deg' = degrees, 'rad' = radians
%
%    Outputs:
%        f      -- true anomaly at epoch
    d2r = pi()/180;

    if(strcmp(angletype,'deg'))
        M = M*d2r;
    end

    M = mod(M,(2*pi()));

    %  Compute eccentric anomaly
    %  A priori estimate

    j = 1;
    E = [];
    if (((-pi() < M)&&(M < 0))||(M > pi()))
        E(end+1) = (M - e);
    else
        E(end+1) = (M + e);
    end

    %  Newton iteration to find eccentric anomaly
    %  Algorithm [goal: find E so f = 0]

    f_E = [E(1) - e*sin(E(1))-M];

    while((abs(f_E(j)) > 1E-13)&&(j <= 100))
        E(end+1) = (E(j) - f_E(j)/(1-e*cos(E(j))));
        j = j + 1;
        f_E(end+1) = (E(j)-e*sin(E(j)) - M);
    end

    % Converged eccentric anomaly [ rad ]    
    E = E(j);

    % True anomaly [ rad ]
    f = mod((2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))),(2*pi()));

    if(strcmp(angletype,'deg'))
        f = f/d2r;
    end
end