function rda = eci2rda(eci)
%    This function computes the spherical equatorial coordinates from the 
%    Earth-centered inertial frame using the Greenwich Mean Standard Time (GMST) and the
%    Earth's eccentricity (ee)
%
%    Parameters:
%        eci -- [ECI position vector [km]]
%
%    Outputs:
%       rda -- [right ascension [ deg ], declination [ deg ], altitude [ km ]]

    % Extract x, y, z from the matrix
    x = eci(:, 1);
    y = eci(:, 2);
    z = eci(:, 3);
    
    % Compute the magnitude of the ECI vectors (distance from Earth's center)
    r = sqrt(x.^2 + y.^2 + z.^2);
    
    % Calculate Declination (Dec) in radians
    Dec = asin(z ./ r);
    
    % Calculate Right Ascension (RA) in radians
    RA = atan2(y, x);
    
    % Ensure RA is in the range [0, 2*pi]
    RA = mod(RA, 2*pi);
    
    % Convert RA from radians to hours
    RA = rad2deg(RA);
    
    % Convert Dec from radians to degrees
    Dec = rad2deg(Dec);
    
    % Calculate Altitude (distance from the Earth's surface)
    earthRadius = 6378.1; % Mean radius of the Earth
    Altitude = r - earthRadius;
    
    % Combine RA, Dec, and Altitude into the output Mx3 matrix
    rda = [RA, Dec, Altitude];
end