function eq_state = ecliptic_to_equatorial(ec_state, obliquity)
% Convert equatorial coordinates to ecliptic coordinates
% 
% Parameters:
%     ec_state  -- Cartesian state in ecliptic frame
%     obliquity -- Obliquity of the central body [ deg ]
% 
% Outputs:
%     eq_state  --  Cartesian position and velocity in equatorial frame, (X, Y, Z, VX, VY, VZ)

    X = ec_state(1);
    Y = ec_state(2);
    Z = ec_state(3);
    VX = ec_state(4);
    VY = ec_state(5);
    VZ = ec_state(6);

    % Convert angles to radians
    obliquity_rad = deg2rad(obliquity);
    
    % Calculate the transformation matrix for the rotation around the x-axis
    R_x = [1, 0, 0;
           0, cos(obliquity_rad), sin(obliquity_rad);
           0, -sin(obliquity_rad), cos(obliquity_rad)];
    
    % Apply the rotation to the position and velocity vectors
    position = [X; Y; Z];
    velocity = [VX; VY; VZ];
    
    eq_position = R_x' * position;
    eq_velocity = R_x' * velocity;

    eq_state = [eq_position; eq_velocity];
end