function [rv_I] = AGI_Qmethod_Inverse_RotateBarycenter_rv_to_J2000_rv(rv_jp, RV_R_barycenter_unitless, const)
% convert state from inertial frame to rotating frame
%   Input:
%             rv_jp: the secondary's position and velocity vector in inertial frame
%             oe: target body's orbita element (a,e,i,Omega,omega) in degree and km
%             const: includes const.primary.mu (the primary gravitational constant), 
%                    const.secondary.mu (the secondary gravitational constant)
%     Output:
%             RV_R_barycenter_unitless: position and velocity vector defined by the 
%             rotating frame of the primary and the secondary with unit length being 
%             the distance between the seondary and the primary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rotating frame unitless to rotating frame with unit
lstar = norm(rv_jp(1:3)); % characteristic length
mustar = const.primary.mu+const.secondary.mu; % G*(characteristic u)
tstar = sqrt(lstar^3/mustar);
R_R_barycenter_unitless = RV_R_barycenter_unitless(1:3); % unitless R in rotating barycenter
V_R_barycenter_unitless = RV_R_barycenter_unitless(4:6); % unitless V in rotating barycenter
R_R_barycenter = R_R_barycenter_unitless*lstar; % unit in km R in rotating barycenter
V_R_barycenter = V_R_barycenter_unitless* lstar/tstar; % unit in km/s V in rotating barycenter
RV_R_barycenter = [R_R_barycenter;V_R_barycenter];

% Take ephemeris of Jupiter position and calculate circular orbit velocity and create DCM
r_jp = rv_jp(1:3); %km
v_jp_ori = rv_jp(4:6); %km

% r3bp requires circle orbit
omega = sqrt(const.primary.mu/ norm(r_jp)^3);
hhat = cross(r_jp,v_jp_ori)/norm(cross(r_jp,v_jp_ori)); % direction of out of plane
omega_vector = omega*hhat;
v_jp_omega = cross(omega_vector,r_jp); % circular assumption velocity vector

% circular assumption to cover v_jp use v_jp_omega from circular
% calculation
v_jp = v_jp_omega; % use circular assumption for velocity vector

% instanteneous rotating axes
xhat = r_jp / norm(r_jp);
zhat = cross(r_jp,v_jp) / norm(cross(r_jp,v_jp));
yhat = cross(zhat,xhat);

% position transfrom related
% rotating frame to inertial frame DCM (right to left)
I_DCM_R = [xhat,yhat,zhat];
% inertial frame to rotating frame DCM (Just the transpose of previous one)
R_DCM_I = I_DCM_R';

% Another way express dxhat dyhat dzhat from AGI page
% benefit no need to mess with angular velocity and its direction
ratio_u = const.secondary.mu / (const.primary.mu+const.secondary.mu);
dxhat = v_jp/norm(r_jp)-r_jp*(r_jp'*v_jp)/norm(r_jp)^3;
dyhat = cross(cross(r_jp,v_jp),v_jp)/(norm(r_jp)*norm(cross(r_jp,v_jp))) -...
    (cross(r_jp,v_jp)/(norm(r_jp)^3*norm(cross(r_jp,v_jp))))'*(cross(cross(r_jp,v_jp),r_jp));
dzhat = [0;0;0];
dQ_matrix = [dxhat,dyhat,dzhat]';
Q = R_DCM_I;
R0 = ratio_u*r_jp;
V0 = ratio_u*v_jp;
QdQQ_matrix = [Q,zeros(3);dQ_matrix,Q];
QdQQ_matrix_inverse = QdQQ_matrix^-1;

% rotating barycenter frame -> shift barycenter to sun-centered inertial rv -> J2000 sun-centered inertial rv
rv_I = QdQQ_matrix_inverse*RV_R_barycenter + [R0;V0];

end

