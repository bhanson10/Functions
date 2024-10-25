function ec_state = equatorial_to_ecliptic(eq_state, obliquity)
% Convert equatorial coordinates to ecliptic coordinates
% 
% Parameters:
%     eq_state  -- Cartesian state in equatorial frame
%     obliquity -- Obliquity of the central body [ deg ]
% 
% Outputs:
%     ec_state  --  Cartesian position and velocity in ecliptic frame, (X, Y, Z, VX, VY, VZ)

    X = eq_state(1);
    Y = eq_state(2);
    Z = eq_state(3);
    VX = eq_state(4);
    VY = eq_state(5);
    VZ = eq_state(6);

    % Convert angles to radians
    obliquity_rad = deg2rad(obliquity);
    
    % Calculate the transformation matrix for the rotation around the x-axis
    R_x = [1, 0, 0;
           0, cos(obliquity_rad), sin(obliquity_rad);
           0, -sin(obliquity_rad), cos(obliquity_rad)];
    
    % Apply the rotation to the position and velocity vectors
    position = [X; Y; Z];
    velocity = [VX; VY; VZ];
    
    ec_position = R_x * position;
    ec_velocity = R_x * velocity;

    ec_state = [ec_position; ec_velocity];
end