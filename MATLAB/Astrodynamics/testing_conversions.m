clear all; close all; clc; format long

%% Initial Conditions 
mu = 398600.4354;
rv = [7007.226732000000; 0; 0; 0; 0.660620156299; 7.550934749459];
rv_cov = diag([(3E-2)^2; (3E-2)^2; (3E-2)^2; (3E-4)^2; (3E-4)^2; (3E-4)^2]);

% mu = 3202.738774922892;
% rv = [1690.689890962806; 0; 0; 0; -0.2954850833646089; -1.352955033637386];
% rv_cov = diag([1^2; 1^2; 1^2; 0.001^2; 0.001^2; 0.001^2]);

%% Converting from RV and RV_cov to COE and COE_cov (working)
coe     = rv2oe(rv,mu,'rad','f','coe')', coe = coe';
coe_cov = rv2oe_cov(rv,rv_cov,mu,'f','coe'), 

%% Converting from RV and RV_cov to EOE and EOE_cov (working)
% eoe     = rv2oe(rv,mu,'rad','f','eoe')', eoe = eoe';
% eoe_cov = rv2oe_cov(rv,rv_cov,mu,[],'eoe'), 

%% Testing Conversion from RV to OE and back (working)
% rv   = rv', rv = rv'; 
% eoe  = rv2oe(rv,mu,'rad','f','eoe')', eoe = eoe';
% rv_r = oe2rv(eoe,mu,'rad','f','eoe')', rv_r = rv_r'; 

%% Testing Conversion from COE to EOE and back (working)
% coe = rv2oe(rv,mu,'rad','f','coe')', coe = coe';
% eoe = oe2oe(coe,'rad','f','coe','eoe')', eoe = eoe';
% coe_r = oe2oe(eoe,'rad','f','eoe','coe')', coe_r = coe_r';

%% Testing Conversion from RV covariance to COE covariance and back (working)
% anomaly = 'M'; % 'f': true anomaly, 'E': eccentric anomaly, 'M': mean anomaly
% 
% rv_cov, 
% coe     = rv2oe(rv,mu,'rad',anomaly,'coe');
% coe_cov = rv2oe_cov(rv,rv_cov,mu,anomaly,'coe'),
% rv_cov  = oe2rv_cov(coe,coe_cov,mu,anomaly,'coe')

%% Testing Conversion from COE covariance to EOE covariance and back (working)
% anomaly = 'f'; % 'f': true anomaly, 'E': eccentric anomaly, 'M': mean anomaly
% 
% coe       = rv2oe(rv,mu,'rad',anomaly,'coe');
% coe_cov   = rv2oe_cov(rv,rv_cov,mu,anomaly,'coe'),
% eoe       = oe2oe(coe,'rad',anomaly,'coe','eoe');
% eoe_cov   = oe2oe_cov(coe,coe_cov,anomaly,'coe','eoe'),
% coe_cov_r = oe2oe_cov(eoe,eoe_cov,anomaly,'eoe','coe'),

%% Testing Conversion from RV covariance to EOE covariance and back (working)
% rv_cov, 
% eoe      = rv2oe(rv,mu,'rad',[],'eoe');
% eoe_cov  = rv2oe_cov(rv,rv_cov,mu,[],'eoe'),
% rv_cov_r = oe2rv_cov(eoe,eoe_cov,mu,[],'eoe'),

