function [RV_R_barycenter_unitless] = AGI_Qmethod_J2000_rv_to_RotateBaryCenter_rv(rv_jp, rv_I, const)
    % J2000 Sun-centered Inertial frame r and v in column vector
    rv_I = rv_I(:);
    % Take ephemeris of Jupiter position and calculate circular orbit velocity and create DCM
    rv_jp = rv_jp(:);
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
    % sun-centered inertial rv -> shift sun-centered inertial rv to barycenter -> rotating barycenter frame
    r_I = rv_I(1:3);
    v_I = rv_I(4:6);
    RV_R_barycenter = QdQQ_matrix*[r_I-R0;v_I-V0];
    % RV_R_barycenter = QdQQ_matrix*[r_I;v_I];
    R_R_barycenter = RV_R_barycenter(1:3);
    V_R_barycenter = RV_R_barycenter(4:6);
    % rotating frame with unit to rotating frame unitless
    lstar = norm(rv_jp(1:3)); % characteristic length
    mustar = const.primary.mu+const.secondary.mu; % G*(characteristic u)
    tstar = sqrt(lstar^3/mustar);
    R_R_barycenter_unitless = R_R_barycenter / lstar;
    V_R_barycenter_unitless = V_R_barycenter * tstar/lstar;
    RV_R_barycenter_unitless = [R_R_barycenter_unitless;V_R_barycenter_unitless];
end