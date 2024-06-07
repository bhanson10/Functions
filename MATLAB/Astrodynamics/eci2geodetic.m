function [long, lat, alt] = eci2geodetic(eci, GMST, angletype)
%    This function computes the geodetic coordinates from the Earth-centered
%    inertial frame using the Greenwich Mean Standard Time (GMST) and the
%    Earth's eccentricity (ee)
%
%    Parameters:
%        eci -- [ECI position vector [km]]
%        GMST -- Greenwich Mean Standard Time [s]
%        angletype --  string indicating desired angletype for inputs
%                      * 'deg' = degrees
%                      * 'rad' = radians
%
%    Outputs:
%        longitude [rad or deg], latitude [rad or deg], altitude [km]

Re = 6371; % Earth radius [km]
ee = 0.017; % Earth eccentricity
if ~exist('angletype', 'var'), angletype = 'rad'; end

Q_eci_ecef = [cos(GMST)  sin(GMST) 0;
              -sin(GMST) cos(GMST) 0;
              0          0         1];
r_ecef = Q_eci_ecef*(eci);
r = sqrt(r_ecef(1)^2+r_ecef(2)^2); 
long = atan2(r_ecef(2), r_ecef(1)); 
lat0 = atan2(r_ecef(3), sqrt(r_ecef(1)^2+r_ecef(2)^2)); diff = inf; 
while diff > 1e-5
    N = Re/sqrt(1-ee^2*sin(lat0)^2); 
    lat1 = atan2(r_ecef(3)+ee^2*N*sin(lat0), sqrt(r_ecef(1)^2+r_ecef(2)^2)); 
    diff = lat1 - lat0; lat0 = lat1; 
end
lat = lat1; 
alt = (sqrt(r_ecef(1)^2+r_ecef(2)^2)/cos(lat)) - N; 

if strcmp(angletype, 'deg')
    long = rad2deg(long); 
    lat = rad2deg(lat); 
end
end