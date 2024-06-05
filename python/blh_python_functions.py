"""
blh_python_functions.py

Author: Benjamin Hanson

Date: 9/7/2023

List of miscellaneous python functions used for UCSD research
"""

import numpy as np
import ssl
ssl._create_default_https_context = ssl._create_unverified_context # This is so we can connect to NASA's JPL Horizons system
from   scipy.stats import norm, t
import spiceypy as cspice
import math
from astropy.time import Time
import scipy
import sys
###########################################################################################################################################
def parse_horizons_ephemeris_koe(file_path):
    """
    Parse Horizons Ephemeris File (Osculating orbital elements)

    Parameters:
        file_path -- path to Horizons ephemeris file 

    Outputs:
        JD  -- Epoch Julian Date,
        TDB -- Barycentric Dynamical Time
        EC  -- Eccentricity
        QR  --  Periapsis distance
        IN  -- Inclination w.r.t. xy-plane (degrees)
        OM  -- Longitude of Ascending Node (degrees)
        W   -- Argument of Perifocus (degrees)
        Tp  -- Periapsis time (user specifies absolute or relative date)
        N   -- Mean motion (degrees/DU)
        MA  -- Mean anomaly (degrees)
        TA  -- True anomaly (degrees)
        A   -- Semi-major axis
        AD  -- Apoapsis distance
        PR  -- Orbital Period) 
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    JD, TDB, EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR = ([] for i in range(14))
    data_start = False

    for line in lines:
        if '$$SOE' in line:
            data_start = True
            line_count = 0
            continue
        if '$$EOE' in line:
            break

        if data_start and line.strip() != '':
            elements = line.strip().split(' ')
            if (line_count % 5 == 0):
                JD.append(float(elements[0]))
                TDB.append(elements[3] + ', ' + elements[4])
            elif(line_count % 5 == 1):
                EC.append(float(elements[1]))
                QR.append(float(elements[3]))
                IN.append(float(elements[5]))
            elif(line_count % 5 == 2):
                OM.append(float(elements[1]))
                W.append(float(elements[4]))
                Tp.append(float(elements[7]))
            elif(line_count % 5 == 3):
                N.append(float(elements[2]))
                MA.append(float(elements[4]))
                TA.append(float(elements[6]))
            else:
                A.append(float(elements[2]))
                AD.append(float(elements[4]))
                PR.append(float(elements[6]))
    
            line_count = line_count + 1 

    return JD, TDB, EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR
###########################################################################################################################################
def parse_horizons_ephemeris_rv(file_path):
    """
    Parse Horizons Ephemeris File (Cartesian states)

    Parameters:
    file_path -- path to Horizons ephemeris file 

    Outputs:
        JD  --  Julian date
        TDB -- Barycentric dynamical time
        RV  -- Carteisan position and velocity {X, Y, Z, VX, VY, VZ} (AU, AU/day) 
    """

    with open(file_path, 'r') as file:
        lines = file.readlines()

    JD, TDB, X, Y, Z, VX, VY, VZ = ([] for i in range(8))
    data_start = False

    for line in lines:
        if '$$SOE' in line:
            data_start = True
            line_count = 0
            continue
        if '$$EOE' in line:
            break

        if data_start and line.strip() != '':
            elements = line.strip().split()
            if (line_count % 3 == 0):
                JD.append(float(elements[0]))
                TDB.append(elements[3] + ' ' + elements[4])
            elif(line_count % 3 == 1):
                X.append(float(elements[2]))
                Y.append(float(elements[3]))
                Z.append(float(elements[4]))
                VX.append(float(elements[5]))
                VY.append(float(elements[6]))
                VZ.append(float(elements[7]))
    
            line_count = line_count + 1 
    RV = np.array([X, Y, Z, VX, VY, VZ]).T

    return JD, TDB, RV

###########################################################################################################################################
def parse_tle_hist_koe(file_path):
    """
    Parse TLE File (provided by Prof. Aaron Rosengren)

    Parameters:
    file_path -- path to Horizons ephemeris file 

    Outputs:
        JD  --  Julian date
        A   -- Semi-major axis [ km ]
        EC  -- Eccentricity [ ]
        IN  -- Inclination w.r.t. ecliptic [ deg ]
        OM  -- Longitude of Ascending Node [ deg ] 
        W   -- Argument of Perifocus [ deg ]
        MA  -- Mean anomaly [ deg ]
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    JD, A, EC, IN, OM, W, MA = ([] for i in range(7))

    for line in lines:
        elements = line.strip().split(',')
        JD.append(float(elements[0]))
        A.append(float(elements[1]))
        EC.append(float(elements[2]))
        IN.append(float(elements[3]))
        OM.append(float(elements[4]))
        W.append(float(elements[5]))
        MA.append(float(elements[6]))

    return JD, A, EC, IN, OM, W, MA
###########################################################################################################################################
def parse_python_rv(file_path, prop_start_JD, prop_end_JD):
    """
    Parse Python-created Ephemeris File (Cartesian states, geocentric, EcJ2000)

    Parameters:
    file_path -- path to Horizons ephemeris file 

    Outputs:
        JD  --  Julian date
        RV  -- Carteisan position and velocity {X, Y, Z, VX, VY, VZ} (KM, KM/sec) 
    """

    with open(file_path, 'r') as file:
        lines = file.readlines()

    JD, X, Y, Z, VX, VY, VZ = ([] for i in range(7))

    for line in lines:
        elements = line.strip().split()
        if((float(elements[0]) >= prop_start_JD)and(float(elements[0]) <= prop_end_JD)):
            JD.append(float(elements[0]))
            X.append(float(elements[1]))
            Y.append(float(elements[2]))
            Z.append(float(elements[3]))
            VX.append(float(elements[4]))
            VY.append(float(elements[5]))
            VZ.append(float(elements[6]))

    RV = np.array([X, Y, Z, VX, VY, VZ]).T

    return JD, RV
###########################################################################################################################################
def parse_GMAT_rv(file_path):
    """
    Parse GMAT-created Ephemeris File (Cartesian states)

    Parameters:
    file_path -- path to Horizons ephemeris file 

    Outputs:
        JD  -- Julian date
        RV  -- Carteisan position and velocity {X, Y, Z, VX, VY, VZ} (KM, KM/sec) 
    """

    with open(file_path, 'r') as file:
        lines = file.readlines()

    JD, X, Y, Z, VX, VY, VZ = ([] for i in range(7))

    count = 0
    for line in lines:
        if count!=0:
            elements = line.strip().split()
            JD.append(float(elements[0])+2430000.0)
            X.append(float(elements[1]))
            Y.append(float(elements[2]))
            Z.append(float(elements[3]))
            VX.append(float(elements[4]))
            VY.append(float(elements[5]))
            VZ.append(float(elements[6]))
        count+=1 

    RV = np.array([X, Y, Z, VX, VY, VZ]).T

    return JD, RV
###########################################################################################################################################
def parse_GMAT_oe(file_path):
    """
    Parse GMAT-created Ephemeris File (Orbital elements)

    Parameters:
    file_path -- path to Horizons ephemeris file 

    Outputs:
        JD   -- Julian date
        A    -- Semi-major axis [ km ]
        EC   -- Eccentricity    [ ]
        INC  -- Inclination relative to the ecliptic [ deg ]
        RAAN -- Inclination relative to the ecliptic [ deg ]
        W    -- Argument of perigee [ deg ]
        TA   -- True anomaly [ deg ]

    """

    with open(file_path, 'r') as file:
        lines = file.readlines()

    JD, A, EC, INC, RAAN, W, TA = ([] for i in range(7))

    count = 0
    for line in lines:
        if count!=0:
            elements = line.strip().split()
            JD.append(float(elements[0])+2430000.0)
            A.append(float(elements[1]))
            EC.append(float(elements[2]))
            INC.append(float(elements[3]))
            RAAN.append(float(elements[4]))
            W.append(float(elements[5]))
            TA.append(float(elements[6]))
        count+=1 

    return JD, A, EC, INC, RAAN, W, TA
###########################################################################################################################################
def equatorial_to_ecliptic(state, obliquity):
    """
    Convert equatorial coordinates to ecliptic coordinates

    Parameters:
        state -- Cartesian state in equatorial frame
        obliquity           -- Obliquity of the central body

    Outputs:
        state  --  Cartesian position and velocity in ecliptic frame, [X, Y, Z, VX, VY, VZ]
    """

    X = state[0]
    Y = state[1]
    Z = state[2]
    VX = state[3]
    VY = state[4]
    VZ = state[5]

    # Convert angles to radians
    obliquity_rad = np.radians(obliquity)
    
    # Calculate the transformation matrix for the rotation around the x-axis
    R_x = np.array([[1, 0, 0],
                    [0, np.cos(obliquity_rad), np.sin(obliquity_rad)],
                    [0, -np.sin(obliquity_rad), np.cos(obliquity_rad)]])
    
    # Apply the rotation to the position and velocity vectors
    position = np.array([X, Y, Z])
    velocity = np.array([VX, VY, VZ])
    
    ecliptic_position = np.dot(R_x, position)
    ecliptic_velocity = np.dot(R_x, velocity)
    
    return np.concatenate((ecliptic_position, ecliptic_velocity),axis=None)

###########################################################################################################################################
def ecliptic_to_equatorial(state, obliquity):
    """
    Convert ecliptic coordinates to equatorial coordinates

    Parameters:
        state -- Cartesian state in equatorial frame
        obliquity           -- Obliquity of the central body

    Outputs:
        state  --  Cartesian position and velocity in equatorial frame, [X, Y, Z, VX, VY, VZ]
    """

    X = state[0]
    Y = state[1]
    Z = state[2]
    VX = state[3]
    VY = state[4]
    VZ = state[5]

    # Convert angles to radians
    obliquity_rad = np.radians(obliquity)
    
    # Calculate the transformation matrix for the rotation around the x-axis
    R_x = np.array([[1, 0, 0],
                    [0, np.cos(obliquity_rad), np.sin(obliquity_rad)],
                    [0, -np.sin(obliquity_rad), np.cos(obliquity_rad)]])
    
    # Apply the rotation to the position and velocity vectors
    position = np.array([X, Y, Z])
    velocity = np.array([VX, VY, VZ])
    
    equatorial_position = np.dot(R_x.T, position)
    equatorial_velocity = np.dot(R_x.T, velocity)

    return np.concatenate((equatorial_position, equatorial_velocity),axis=None)

###########################################################################################################################################
def rv2koe(state, mu, angletype):
    """
    This function computes the classical (Keplerian) orbit
    elements given the Cartesian position and velocity vectors, the
    gravitational parameter of the attracting body, and the desired angletype

    Parameters:
        r -- Cartesian position vector [ km ]
        v -- Cartesian velocity vector [ km/s ]
        mu -- gravitational parameter of attracting body [ km^3/s^2 ]
        angletype --  string indicating desired angletype for inputs
                      'deg' = degrees, 'rad' = radians

    Outputs:
        koe -- vector which contains the classical orbit elements
               [ a; e; i; Om; w; M ]
    """
    r = state[0:3]
    v = state[3:6]

    koe = np.zeros(6)
    r2d = 180/math.pi 
    I = [1, 0, 0]
    J = [0, 1, 0]
    K = [0, 0, 1]

    rhat = r/np.linalg.norm(r) # position unit vector [km]
    h    = np.cross(r, v)      # angular momentum vector
    hhat = h/np.linalg.norm(h) # normalized angular momentum 
    nhat = np.cross(K,h)/np.linalg.norm(np.cross(K,h)) #  normalized ascending node vector

    # Eccentricity
    e = (1/mu)*np.cross(v,h) - rhat # Eccentricity vector
    koe[1] = np.linalg.norm(e)

    energy = (1/2)*np.dot(v,v) - mu/np.linalg.norm(r) # energy km^2/s^2
    # If energy < 0, the orbit is closed (periodic)

    # Semi-major axis (a) and parameter (p)

    if(koe[1] != 1):
        koe[0] = -mu/(2*energy)
        p      = koe[0]*(1-koe[1]**2)
    else:
        koe[0] = np.inf 
        p      = np.linalg.norm(h)**2/mu

    # Inclination (i) of orbit
    koe[2] = np.arccos(np.dot(K, hhat))
    # If i < 90 deg, the elliptical orbit is a direct (prograde) orbit

    # Right ascension of the ascending node (Omega)
    koe[3] = np.arctan2(np.dot(J,nhat),np.dot(I,nhat))%(2*math.pi)

    # Argument of periapsis (w)
    koe[4] = np.arctan2(np.dot(hhat, np.cross(nhat,e)), np.dot(nhat,e))%(2*math.pi)

    # True anomaly (f) at epoch [rad]
    f = np.arctan2(np.dot(hhat, np.cross(e,r)),np.dot(e,r))%(2*math.pi)

    # Eccentric anomaly (E) at epoch [rad]
    E = 2*np.arctan2(np.sqrt(1-koe[1])*np.tan(f/2),np.sqrt(1+koe[1]))

    # Mean anomaly (M) at epoch [rad]
    koe[5] = (E-koe[1]*np.sin(E))%(2*math.pi)

    if(angletype == 'deg'):
        koe[2:7]*=r2d

    return koe

###########################################################################################################################################
def oe2rv(oe, mu, oetype, angletype):
    """
    This function computes the orbit elements 
    given the Cartesian position and velocity vectors, the gravitational
    parameter of the attracting body, and the desired angletype

    Parameters:
        oe        -- vector which contains the nonsingular orbit elements
                        [ a; th; in; q1; q2; Om ]
                            a  - semimajor axis [ km ]
                            th - argument of latitude [ rad ]
                            in - inclination [ rad ]
                            q1 - ec*cos(w) [ ]
                            q2 - ec*sin(w) [ ]
                            Om - longitude of the ascending node [ rad ]
                    OR the classical orbit elements
                        [ a; ec; in; Om; w; M ]
                            ec - eccentricity [ ]
                            w  - argument of periapsis [ rad ]
                            M0 - mean anomaly at epoch [ rad ]
        mu        -- gravitational parameter of attracting body [ km^3/s^2 ]
        oetype    -- string indicating orbit element input type
                        'noe' = nonsingular, 'coe' = classical (Keplerian)
        angletype -- string indicating desired angletype for inputs
                        'deg' = degrees, 'rad' = radians

    Outputs:
        rv     -- (Cartesian position vector [ km ], Cartesian velocity vector [ km/s ])  
    """

    d2r = math.pi/180

    if(oetype == 'noe'):

        # Nonsingular orbit elements
        a   = oe[0] # semi-major axis [ km ]
        th  = oe[1] # argument of latitude [ rad ]
        inc = oe[2] # inclination [ rad ]
        q1  = oe[3] # ec*cos(w) [ ]
        q2  = oe[4] # ec*sin(w) [ ]
        Om  = oe[5] # long of the ascending node [ rad ]

        if(angletype == 'deg'):
            th*=d2r
            inc*=d2r
            Om*=d2r
        
        # Mod out by 2pi
        th%=(2*math.pi)
        inc%=math.pi
        Om%=(2*math.pi)

        # Extract classical orbit elements 
        ec = np.sqrt(q1**2 + q2**2) # eccentricity [ ]
        w  = np.arctan2(q2, q1)     # argument of periapsis [ rad ]

        # True anomaly and eccentric anomaly at epoch [ rad ]
        f = (th - w)%(2*math.pi)
        E = (2*np.arctan(np.sqrt((1 - ec)/(1 + ec))*np.tan(f/2)))%(2*math.pi)

        # Mean anomaly at epoch [ deg ]
        M = E - ec*np.sin(E)

    else:

        #  Classical orbit elements
        a   = oe[0]    
        ec  = oe[1]   
        inc = oe[2]   
        Om  = oe[3]   
        w   = oe[4]  
        M   = oe[5]   
            
        if(angletype == 'deg'):
            inc*=d2r
            Om*=d2r
            w*=d2r
            M*=d2r
        
        #  Mod out by pi or 2pi
        inc%=math.pi
        Om%=(2*math.pi) 
        w%=(2*math.pi) 
        M%=(2*math.pi)  

    if (a > 0):  # ----- eccentric orbit ----- #
        # True anaomaly [ rad ]
        f = KepEqn(M, ec, 'rad')
    else:        # ----- hyperbolic orbit ----- #
        
        #  Compute hyperbolic anomaly
        #  A priori estimate
        j = 0
        H = [M]

        #  Newton iteration to find hyperbolic anomaly
        #  Algorithm [goal: find H so f = 0]
        f_H = [ec*np.sinh(H[j]) - H[j] - M]
        while abs(f_H[j]) > 1e-13:
            H.append(H[j] - f_H[j]/(ec*np.cosh(H[j]) - 1))
            j+=1
            f_H.append(ec*np.sinh(H[j]) - H[j] - M)

        #  Converged eccentric anomaly [ rad ]
        H = H[j]

        #  True anomaly [ rad ]
        f = (2*np.arctan(np.sqrt((ec + 1)/(ec - 1))*np.tanh(H/2)))%(2*math.pi)

    # Argument of latitude [ rad ]
    th = w + f 

    # Magnitude of position vector [ km ]
    r = a*(1-ec**2)/(1+ec*np.cos(f)) # trajectory equation

    # Magnitude of velocity vector [ km/s ]
    v = np.sqrt(mu*(2/r-1/a)) # vis-viva equation

    # Ascending node vector
    nhat = np.array([np.cos(Om), np.sin(Om), 0])
    rT   = np.array([-np.cos(inc)*np.sin(Om), np.cos(inc)*np.cos(Om), np.sin(inc)])

    gamma = np.arctan2(ec*np.sin(f), 1+ec*np.cos(f)) # [ rad ]

    # Normalized position and velocity vectors
    rhat = np.cos(th)*nhat + np.sin(th)*rT # [ km ]
    vhat = np.sin(gamma - th)*nhat + np.cos(gamma - th)*rT # [ km/s ]

    # Position and velocity vectors
    r = r*rhat # [ km ]
    v = v*vhat # [ km/s ]

    rv = np.concatenate((r,v), axis=None)

    return rv
###########################################################################################################################################
def newtonRaphson(EA, M, e):
    EPSILON = 0.001
    h = true_anomaly(EA, M, e)/der_true_anomaly(EA, e);
    while (abs(h) >= EPSILON):
        h = true_anomaly(EA, M, e)/der_true_anomaly(EA, e);
        EA -= h
    return EA 
###########################################################################################################################################
def true_anomaly(EA, M, e):
    f = EA - e*np.sin(EA) - M
    return f
###########################################################################################################################################
def der_true_anomaly(EA, e):
    df = 1 - e*np.cos(EA)
    return df
###########################################################################################################################################
def KepEqn(M, e, angletype):
    """
    Solves Kepler's equation 

    Parameters:
    M          -- mean anomaly at epoch 
    e          -- eccentricity
    angletype  -- string indicating desired angletype for inputs
                    'deg' = degrees, 'rad' = radians

    Outputs:
        f      -- true anomaly at epoch
    """

    d2r = math.pi/180

    if(angletype=='deg'):
        M = M*d2r

    M%=(2*math.pi)

    #  Compute eccentric anomaly
    #  A priori estimate

    j = 0
    E = []
    if (((-math.pi < M)and(M < 0))or(M > math.pi)):
        E.append(M - e)
    else:
        E.append(M + e)

    #  Newton iteration to find eccentric anomaly
    #  Algorithm [goal: find E so f = 0]

    f_E = [E[0] - e*np.sin(E[0])-M]

    while((abs(f_E[j]) > 1E-13)and(j <= 100)):
        E.append(E[j] - f_E[j]/(1-e*np.cos(E[j])))
        j+=1
        f_E.append(E[j]-e*np.sin(E[j]) - M)

    # Converged eccentric anomaly [ rad ]    
    E = E[j]

    # True anomaly [ rad ]
    f = (2*np.arctan(np.sqrt((1 + e)/(1 - e))*np.tan(E/2)))%(2*math.pi)

    if(angletype=='deg'):
        f/=d2r

    return f
###########################################################################################################################################
def convert_time(t, out_format, n):
    """
    This function converts time to different formats using SPICE

    Parameters:
        t          -- time, in whatever format (UTC, etc.) 
        output_str -- format of output time (J, C)
        n          -- number of decimals out

    Outputs:
        t_out      -- converted time in output_str format
    """

    if((type(t)== int)or(type(t)==float)):
        time_str = "JD" + str(t)
    else:
        time_str = t
    
    try:
        et = cspice.str2et(time_str)
    except:
        print("Invalid input time format.")
        return
    
    try:
        t_out = cspice.et2utc(et, out_format, n)
    except:
        print("Invalid input time format.")
        return
    
    if (out_format == 'J'):
        t_out = float(str.split(t_out)[1])

    return t_out
###########################################################################################################################################
class const:
    def __init__(self):
        self.primary_mu   = 0 
        self.secondary_mu = 0 
        self.mu           = 0
        self.LU           = 0
        self.TU           = 0
        self.per_of_LHR   = 0
        self.hill_r       = 0

###########################################################################################################################################
class trajectory:
    def __init__(self):
        self.JD       = []
        self.geo_RV   = [] 
        self.syn_RV   = []
        self.resid_geo_RV = []
        self.pos_geo_err  = []
        self.vel_geo_err  = []
        self.resid_syn_RV = []
        self.pos_syn_err  = []
        self.vel_syn_err  = []
###########################################################################################################################################
def AGI_Qmethod_J2000_rv_to_RotateBaryCenter_rv(rv_sec, rv_I, const):
    r_sec = rv_sec[0:3]
    v_sec_ori = rv_sec[3:6]
    
    omega = np.sqrt(const.primary_mu / np.linalg.norm(r_sec)**3)
    hhat = np.cross(r_sec, v_sec_ori) / np.linalg.norm(np.cross(r_sec, v_sec_ori))
    omega_vector = omega * hhat
    v_sec_omega = np.cross(omega_vector, r_sec)
    
    v_sec = v_sec_omega
    
    xhat = r_sec / np.linalg.norm(r_sec)
    zhat = np.cross(r_sec, v_sec) / np.linalg.norm(np.cross(r_sec, v_sec))
    yhat = np.cross(zhat, xhat)
    
    I_DCM_R = np.column_stack((xhat, yhat, zhat))
    R_DCM_I = I_DCM_R.T
    
    ratio_u = const.secondary_mu / (const.primary_mu + const.secondary_mu)
    dxhat = v_sec / np.linalg.norm(r_sec) - r_sec * np.dot(np.array(r_sec).T, v_sec) / np.linalg.norm(r_sec)**3
    dyhat = np.cross(np.cross(r_sec, v_sec), v_sec) / (np.linalg.norm(r_sec) * np.linalg.norm(np.cross(r_sec, v_sec))) - np.dot(np.cross(r_sec, v_sec).T / (np.linalg.norm(r_sec)**3 * np.linalg.norm(np.cross(r_sec, v_sec))),
             np.cross(np.cross(r_sec, v_sec), r_sec))
    dzhat = np.array([0, 0, 0])
    dQ_matrix = np.column_stack((dxhat, dyhat, dzhat)).T
    
    Q = R_DCM_I
    R0 = ratio_u * r_sec
    V0 = ratio_u * v_sec
    QdQQ_matrix = np.column_stack((Q, np.zeros((3, 3))))
    QdQQ_matrix = np.vstack((QdQQ_matrix, np.column_stack((dQ_matrix, Q))))
    r_I = rv_I[0:3]
    v_I = rv_I[3:6]

    RV_R_barycenter = np.dot(QdQQ_matrix, np.concatenate((r_I - R0, v_I - V0)))
    R_R_barycenter = RV_R_barycenter[0:3]
    V_R_barycenter = RV_R_barycenter[3:6]
    
    lstar = np.linalg.norm(rv_sec[0:3])
    mustar = const.primary_mu + const.secondary_mu
    tstar = np.sqrt(lstar**3 / mustar)
    
    R_R_barycenter_unitless = R_R_barycenter / lstar
    V_R_barycenter_unitless = V_R_barycenter * tstar / lstar
    
    return np.concatenate([R_R_barycenter_unitless, V_R_barycenter_unitless])

###########################################################################################################################################
def geo_2_syn(JD, geo_RV, c):
    """
    Convert geocentric coordinates to synodic coordinates 

    Parameters:
        JD     -- Julian date 
        geo_RV -- Geocentric position/velocity coordinates, [n, 6] matrix

    Outputs:
        syn_RV -- Synodic position/velocity coordinates, [n, 6] matrix
    """

    cspice.furnsh('/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Kernels/planets/de441_part-1.bsp')
    cspice.furnsh('/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Kernels/planets/de441_part-2.bsp')
    cspice.furnsh('/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Kernels/planets/Kernels/naif0012.tls')

    syn_RV = np.zeros((len(JD),6))
    if (len(JD) == 1):
        rv_I = np.array(geo_RV)
        rv_sec = np.array(cspice.spkezr('Moon', cspice.str2et('JD'+str(JD[j])), 'ECLIPJ2000', 'none', 'Earth')[0])
        syn_RV[0] = AGI_Qmethod_J2000_rv_to_RotateBaryCenter_rv(rv_sec, rv_I, c)
    else:
        for j in range(len(JD)):
                rv_I = np.array([geo_RV[j][0], geo_RV[j][1], geo_RV[j][2], geo_RV[j][3], geo_RV[j][4], geo_RV[j][5]])
                rv_sec = np.array(cspice.spkezr('Moon', cspice.str2et('JD'+str(JD[j])), 'ECLIPJ2000', 'none', 'Earth')[0])
                syn_RV[j] = AGI_Qmethod_J2000_rv_to_RotateBaryCenter_rv(rv_sec, rv_I, c)

    return syn_RV

###########################################################################################################################################
def AGI_Qmethod_Inverse_RotateBarycenter_rv_to_J2000_rv(rv_sec, RV_R_barycenter_unitless, const):

    # rotating frame unitless to rotating frame with unit
    lstar = np.linalg.norm(rv_sec[0:3]) # characteristic length
    mustar = const.primary_mu + const.secondary_mu # G*(characteristic u)
    tstar = np.sqrt(lstar**3 / mustar) # characteristic time
    R_R_barycenter_unitless = RV_R_barycenter_unitless[0:3]  # unitless R in rotating barycenter
    V_R_barycenter_unitless = RV_R_barycenter_unitless[3:6]  # unitless V in rotating barycenter
    R_R_barycenter = R_R_barycenter_unitless * lstar         # unit in km R in rotating barycenter
    V_R_barycenter = V_R_barycenter_unitless * lstar / tstar # unit in km/s V in rotating barycenter
    RV_R_barycenter = np.concatenate([R_R_barycenter,V_R_barycenter])

    # Take ephemeris of Jupiter position and calculate circular orbit velocity and create DCM
    r_sec = rv_sec[0:3]     # km
    v_sec_ori = rv_sec[3:6] # km

    # r3bp requires circle orbit
    omega = np.sqrt(const.primary_mu / np.linalg.norm(r_sec)**3)
    hhat = np.cross(r_sec, v_sec_ori) / np.linalg.norm(np.cross(r_sec, v_sec_ori)) # direction of out of plane
    omega_vector = omega * hhat
    v_sec_omega = np.cross(omega_vector, r_sec) # circular assumption velocity vector

    # circular assumption to cover v_sec use v_sec_omega from circular calculation
    v_sec = v_sec_omega # use circular assumption for velocity vector

    # instanteneous rotating axes
    xhat = r_sec / np.linalg.norm(r_sec)
    zhat = np.cross(r_sec, v_sec) / np.linalg.norm(np.cross(r_sec, v_sec))
    yhat = np.cross(zhat, xhat)

    # position transfrom related rotating frame to inertial frame DCM (right to left)
    I_DCM_R = np.column_stack((xhat, yhat, zhat))
    # inertial frame to rotating frame DCM (Just the transpose of previous one)
    R_DCM_I = I_DCM_R.T

    # Another way express dxhat dyhat dzhat from AGI page benefit no need to mess with angular velocity and its direction
    ratio_u = const.secondary_mu / (const.primary_mu + const.secondary_mu)
    dxhat = v_sec / np.linalg.norm(r_sec) - r_sec * np.dot(np.array(r_sec).T, v_sec) / np.linalg.norm(r_sec)**3
    dyhat = np.cross(np.cross(r_sec, v_sec), v_sec) / (np.linalg.norm(r_sec) * np.linalg.norm(np.cross(r_sec, v_sec))) - np.dot(np.cross(r_sec, v_sec).T / (np.linalg.norm(r_sec)**3 * np.linalg.norm(np.cross(r_sec, v_sec))),
             np.cross(np.cross(r_sec, v_sec), r_sec))
    dzhat = np.array([0, 0, 0])
    dQ_matrix = np.column_stack((dxhat, dyhat, dzhat)).T
    Q = R_DCM_I
    R0 = ratio_u * r_sec
    V0 = ratio_u * v_sec
    QdQQ_matrix = np.column_stack((Q, np.zeros((3, 3))))
    QdQQ_matrix = np.vstack((QdQQ_matrix, np.column_stack((dQ_matrix, Q))))
    QdQQ_matrix_inverse = np.linalg.inv(QdQQ_matrix)

    # rotating barycenter frame -> shift barycenter to sun-centered inertial rv -> J2000 sun-centered inertial rv
    rv_I = np.dot(QdQQ_matrix_inverse,RV_R_barycenter) + np.concatenate([R0,V0])

    return rv_I
###########################################################################################################################################
def syn_2_geo(JD, syn_RV, c):
    """
    Convert geocentric coordinates to synodic coordinates 

    Parameters:
        JD     -- Julian date 
        syn_RV -- Synodic position/velocity coordinates, [n, 6] matrix

    Outputs:
        geo_RV -- Geocentric position/velocity coordinates, [n, 6] matrix
    """

    cspice.furnsh('/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Kernels/planets/de441_part-1.bsp')
    cspice.furnsh('/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Kernels/planets/de441_part-2.bsp')
    cspice.furnsh('/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Kernels/leap-seconds/naif0012.tls')

    syn_RV = np.zeros((len(JD),6))
    for j in range(len(JD)):
            RV_R_barycenter_unitless = np.array([syn_RV[j][0], syn_RV[j][1], syn_RV[j][2], syn_RV[j][3], syn_RV[j][4], syn_RV[j][5]])
            rv_sec = np.array(cspice.spkezr('Moon', cspice.str2et('JD'+str(JD[j])), 'ECLIPJ2000', 'none', 'Earth')[0])
            syn_RV[j] = AGI_Qmethod_Inverse_RotateBarycenter_rv_to_J2000_rv(rv_sec, RV_R_barycenter_unitless, const)

    return syn_RV

###########################################################################################################################################
def window_data(JD, data, start_JD, end_JD, resid_step):
    """
    Constrict data set to an observation window

    Parameters:
        JD           -- Julian date data
        data         -- Data corresponding to Julian date length
        start_JD     -- Start of observation window (JD)
        end_JD       -- End of observation window (JD)
        resid_step   -- Time between measurements (hours)

    Outputs:
        window_data  -- Data corresponding to observation window
    """


    window_data = []

    count = 0
    for i in range(len(JD)):
        if((JD[i] >= start_JD)and(JD[i] <= end_JD)):
            window_data.append(data[i])
            count += 1

    return np.array(window_data)

###########################################################################################################################################
def get_orbit_properties(JD, RV, c):
    """
    Constrict data set to an observation window

    Parameters:
        JD -- Julian date data
        RV -- Cartesian position/velocity Data corresponding to Julian date length
        c  -- System constants 

    Outputs:
        dist_to_moon -- Distance to moon in synodic frame
        perilune     -- Times of perilune (JD)
        hill_enters  -- Times entering LHR defined by c.hill_r (JD)
        hill_exits   -- Times exiting LHR defined by c.hill_r (JD)
    """


    dist_to_moon = np.zeros(len(JD))
    for i in range(len(JD)):
        dist_to_moon[i] = ((RV[i][0]-(1-c.mu)*c.LU)**2+(RV[i][1])**2+(RV[i][2])**2)**(1/2)

    perilune_idx = np.where((dist_to_moon[:-2] > dist_to_moon[1:-1]) & (dist_to_moon[2:] > dist_to_moon[1:-1]))[0] + 1

    perilune = []
    for i in perilune_idx:
        perilune.append(JD[i])

    hill_enters = []
    hill_exits = []

    for i in range(len(dist_to_moon)):
        if(i==0):
            if(dist_to_moon[0] <= c.hill_r):
                hill_enters.append(JD[i])
                inside = True
            else:
                inside = False
        elif(i==len(dist_to_moon)-1):
            if(len(hill_enters) > len(hill_exits)):
                hill_exits.append(JD[i])
        else:
            if(dist_to_moon[i] > c.hill_r):
                if(inside==True):
                    hill_exits.append(JD[i])
                    inside = False
            else:
                if(inside==False):
                    hill_enters.append(JD[i])
                    inside = True

    return dist_to_moon, perilune, hill_enters, hill_exits

###########################################################################################################################################
def get_assist_trajectories(start, end, step, ROOT_PATH, e_w_JD, e_w_syn_RV, e_w_geo_RV, e_TDB, c):
    """
    Constrict data set to an observation window

    Parameters:
        start      -- Start IC of ASSIST trajectory in observation window
        end        -- End IC of ASSIST trajectory in observation window
        step       -- Steps between IC files
        ROOT_PATH  -- Path where ASSIST propagated files are 
        e_w_JD     -- Ephemeris observation window time set (JD)
        e_w_syn_RV -- Ephemeris synodic Cartesian position/velocity in observation window
        e_w_syn_RV -- Ephemeris geocentric Cartesian position/velocity in observation window
        e_TDB      -- Barydynamical time
        c          -- const()

    Outputs:
        assist_trajectories -- List of trajectories() with ICs propagated by ASSIST
    """

    IC_list = [start + i * step for i in range((end - start) // step + 1)]

    assist_trajectories = []
    for i in IC_list:
        traj = trajectory()
        traj.JD, traj.geo_RV = parse_python_rv(ROOT_PATH + '/ASSIST Prop/ASSIST_' + e_TDB[i] + '.txt', e_w_JD[0], e_w_JD[-1])
        traj.syn_RV = geo_2_syn(traj.JD, traj.geo_RV, c)

        traj.resid_geo_RV, traj.resid_syn_RV = (np.zeros((len(traj.JD),6)) for k in range(2))
        traj.pos_geo_err, traj.vel_geo_err, traj.pos_syn_err, traj.vel_syn_err = (np.zeros((len(traj.JD))) for k in range(4))
        for j in range(len(traj.JD)):
            if j==0:
                e_start_idx = np.where(e_w_JD == traj.JD[0])[0]
            else:
                e_start_idx += 1
            
            traj.resid_geo_RV[j] = abs(traj.geo_RV[j] - e_w_geo_RV[e_start_idx])
            traj.pos_geo_err[j]  = np.linalg.norm(traj.resid_geo_RV[j][0:3])
            traj.vel_geo_err[j]  = np.linalg.norm(traj.resid_geo_RV[j][3:6])

            traj.resid_syn_RV[j] = abs(traj.syn_RV[j] - e_w_syn_RV[e_start_idx])
            traj.pos_syn_err[j]  = np.linalg.norm(traj.resid_syn_RV[j][0:3])
            traj.vel_syn_err[j]  = np.linalg.norm(traj.resid_syn_RV[j][3:6])

        assist_trajectories.append(traj)

    return assist_trajectories
    
###########################################################################################################################################
def get_maneuvers(trajectories, hill_enters, hill_exits, vel_TOL):
    """
    Constrict data set to an observation window

    Parameters:
        trajectories -- Propagated trajectories() 
        hill_enters  -- List of times marking entry into hill sphere by S/C (JD)
        hill_exits   -- List of times marking exit of hill sphere by S/C (JD)
        vel_TOL      -- Tolerance of detected velocity residual jump (km/s)
        
    Outputs:
        assist_trajectories -- List of trajectories() with ICs propagated by ASSIST
    """

    maneuvers = []
    maneuver_count = []

    for i in range(len(trajectories)):
        traj = trajectories[i]
        for j in range(len(traj.vel_syn_err)-1):
            step = traj.vel_syn_err[j+1]-traj.vel_syn_err[j]
            if(step > vel_TOL):
                maneuver_time = traj.JD[j]
                in_LHR = False
                for k in range(len(hill_enters)):
                    if((maneuver_time >= hill_enters[k])and(maneuver_time <= hill_exits[k])):
                        in_LHR=True
                if(not(in_LHR)):        
                    if(traj.JD[j] not in maneuvers):
                        maneuvers.append(maneuver_time)
                        maneuver_count.append(1)
                    else:
                        idx = maneuvers.index(traj.JD[j])
                        maneuver_count[idx] += 1
                    break
    return maneuvers, maneuver_count

###########################################################################################################################################
def utc_2_tdb(UTCJD):
    """
    Constrict data set to an observation window

    Parameters:
        UTCJD -- Julian Date in UTC time zone
        
    Outputs:
        TDBJD -- Julian Date in TDB time zone
    """

    utc_time = Time(UTCJD, format='jd', scale='utc')
    tdb_time = utc_time.tdb
    TDBJD = tdb_time.jd

    return TDBJD
###########################################################################################################################################
def tdb_2_utc(TDBJD):
    """
    Constrict data set to an observation window

    Parameters:
        TDBJD -- Julian Date in TDB time zone
        
    Outputs:
        UTCJD -- Julian Date in UTC time zone
    """

    tdb_time = Time(TDBJD, format='jd', scale='tdb')
    utc_time = tdb_time.utc
    UTCJD = utc_time.jd

    return UTCJD
###########################################################################################################################################
def convert_tle_file_utc_2_tdb(filepath):
    """
    Convert the epoch format of TLE file from UTC to TDB

    Parameters:
        filepath -- Path to file that is being edited
        
    """

    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Create a list to store modified lines
    modified_lines = []

    # Loop through the lines and edit the UTC JD values
    for line in lines:
        data = line.strip().split(',')
        TDBJD = utc_2_tdb(float(data[0])+ 2400000.5)
        data[0] = str(TDBJD)
        modified_line = ','.join(data) + '\n'
        modified_lines.append(modified_line)

    # Open the file in write mode and write the modified lines back
    with open(filepath, 'w') as file:
        file.writelines(modified_lines)

    return 
###########################################################################################################################################
def EKF_D(x0, xest0, Pest0, dt, T, Q, R, f, h, F=[], H=[], const=[]):

    # Perform Discrete Extended Kalman Filter estimation given dynamics and measurement models

    # Parameters:
    #     x0:      Initial true state
    #     xest0:   Initial estimate
    #     Pest0:   Initial covariance
    #     T:       Period of propagation
    #     dt:      Time step or [True timestep, Measurement timestep]
    #     Q:       Process noise matrix
    #     R:       Measurement noise matrix
    #     f:       Dynamics model
    #     h:       Measurement model
    #     F:       Dynamics Jacobian (optional, finite difference approximation if not included)
    #     H:       Measurement Jacobian (optional, finite difference approximation if not included)
    #     const:   Miscellaneous constants (optional)

    # Outputs:
    #     x:       True states
    #     xest:    Estimated states
    #     z:       Measurements
    #     Pest:    Estimated covariances
    #     tspan_x: True time span
    #     tspan_z: Measurement time span

    # Checks and Balances
    [rows, columns] = np.shape(Q)
    if(rows!=columns):
        print("ERROR: Process noise matrix is not square.")
        sys.exit()

    if(len(xest0)!=len(x0)):
        print("ERROR: Initial true state and initial estimate state dimensions do not match.")
        sys.exit()

    if(len(xest0)!=rows):
        print("ERROR: State vector and process noise matrix are not compatible dimensions.")
        sys.exit()

    [rows, columns] = np.shape(R)
    if(rows!=columns):
        print("ERROR: Measurement noise matrix is not square.")
        sys.exit()

    [rows, columns] = np.shape(Pest0)
    if(rows!=columns):
        print("ERROR: Covariance matrix is not square.")
        sys.exit()

    if(len(xest0)!=rows):
        print("ERROR: State vector and covariance matrix are not compatible dimensions.")
        sys.exit()

    if not(np.allclose(Pest0, Pest0.T, rtol=1e-05, atol=1e-05)):
        print("ERROR: Covariance matrixx is not symmetric.")
        sys.exit()

    if not(callable(f)):
        print("ERROR: Input f must be callable.")
        sys.exit()

    if not(callable(h)):
        print("ERROR: Input h must be callable.")
        sys.exit()

    # Timespans
    if(type(dt) is list)or(type(dt) is np.ndarray):
        if(float(dt[1])/dt[0] % 1 != 0):
            print("Measurement time step is not a multiple of true time step. Rounding to nearest multiple.")
            dt[1] = round(dt[1]/dt[0])*dt[0]; 
            
        tspan_x = np.arange(0,T+dt[0],dt[0])
        tspan_z = np.arange(dt[1],T+dt[1],dt[1])
        dt_ratio = dt[1]/dt[0]
        dt = dt[0]
    else:
        tspan_x = np.arange(0,T+dt,dt)  # Defining true state timespan
        tspan_z = np.arange(dt,T+dt,dt) # Defining measurement timespan
        dt_ratio = 1

    # Approximating Jacobians 
    if(F == []):
        def F(x, dt, const):
            dx = 1e-8
            n = len(x)
            func = f(x, dt, const)
            jac = np.zeros((n, n))
            for j in range(n):  
                Dxj = (abs(x[j])*dx if x[j] != 0 else dx)
                x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
                jac[:, j] = (f(x_plus, dt, const) - func)/Dxj
            return jac

    if(H == []):
        def H(x):
            dx=1e-8
            n = len(x)
            func = h(x)
            jac = np.zeros((n, n))
            for j in range(n):  
                Dxj = (abs(x[j])*dx if x[j] != 0 else dx)
                x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
                jac[:, j] = (h(x_plus) - func)/Dxj
            return jac

    d1      = len(Q)                      # Dimension of process noise matrix
    d2      = len(R)                      # Dimension of measurement noise matrix
    n1      = len(tspan_x)                # Length of true timespan
    n2      = len(tspan_z)                # Length of measurement timespan
    x       = np.full((n1,d1),np.nan)     # True states 
    x[0]    = x0                          # Initializing first true state
    xest    = np.full((n1,d1),np.nan)     # Estimated states
    xest[0] = xest0                       # Initialized first estimated state
    Pest    = np.full((n1,d1,d1),np.nan)  # Estimated covariance
    Pest[0] = Pest0                       # Initializing first estimated covariance
    z       = np.full((n2,d2),np.nan)     # Measurements
    count   = 0
    for i in range(1,n1):

        # Initialization
        w    = np.linalg.multi_dot([np.sqrt(Q),np.random.randn(d1)])    # Process noise
        x[i] = f(x[i-1],dt,const) + w                                   # True state

        # Prediction
        xpred = f(xest[i-1],dt,const)                                                                             # Predicted state
        Ppred = np.linalg.multi_dot([F(xest[i-1],dt,const),Pest[i-1],np.transpose(F(xest[i-1],dt,const))]) + Q    # Predicted covariance
        zpred = h(xpred)                                                                                          # Predicted measurement

        # Correction
        if((i) % dt_ratio == 0):
            v       = np.matmul(np.sqrt(R), np.random.randn(d2))                                                                   # Measurement noise
            z[count]  = h(x[i]) + v                                                                                                # Actual measurement
            K       = np.matmul(Ppred,np.transpose(H(xpred)),np.linalg.inv(np.matmul(H(xpred),Ppred,np.transpose(H(xpred)))+R))    # Kalman Gain
            xest[i] = xpred + np.matmul(K,z[count]-zpred)                                                                          # Estimated state
            Pest[i] = np.matmul((np.identity(d1)-np.matmul(K,H(xpred))),Ppred)                                                     # Estimated covariance
            count += 1
        else:
            xest[i] = xpred                                                                
            Pest[i] = Ppred       

    return x, xest, z, Pest, tspan_x, tspan_z
###########################################################################################################################################
def EKF(xest0, Pest0, dt, z, Q, R, f, time, h = [], F=[], H=[], method = 'RK4', const=[]):

    # Perform Extended Kalman Filter estimation given dynamics and measurement models

    # Parameters:
    #     xest0:   Initial estimate (required)
    #     Pest0:   Initial covariance (required)
    #     dt:      Initial time step (required)
    #     z:       Measurements with epochs (required)
    #       - Measurements can be in the following two formats
    #           * z = [T] where T is the period, meaning no measurements
    #           * z = [[t1, t2, t3, ..., T],[z1,z2,z3,...,zn]] where ti are the measurement epochs and zi are the measurements
    #     Q:       Process noise matrix (required)
    #     R:       Measurement noise matrix (required)
    #     f:       Dynamics model (continuous-time, required, function handle)
    #     time:    Time frame of dynamics model, either "CT" or "DT" (required)
    #     h:       Measurement model (optional if no measurements, function handle)
    #     F:       Dynamics Jacobian (continuous-time, optional, function handle)
    #     H:       Measurement Jacobian (optional, function handle)
    #     method:  Time-marching method (optional)
    #       - 'EE':    Explicit Euler - Dynamics Jacobian is used
    #       - 'RK4':   Runge-Kutta 4 (default) - Dynamics Jacobian is estimated
    #       - 'RK45':  Adaptive Runge-Kutta 4/5 - Dynamics Jacobian is estimated
    #     const:   Miscellaneous constants (optional)

    # Outputs:
    #     xest:    Estimated states
    #     Pest:    Estimated covariances
    #     tspan:   Time span
    
    # Note: RK45 does not work currently
        
    # Checks and Balances
    [rows, columns] = np.shape(Q)
    if(rows!=columns):
        print("ERROR: Process noise matrix is not square.")
        sys.exit()

    if(len(xest0)!=rows):
        print("ERROR: State vector and process noise matrix are not compatible dimensions.")
        sys.exit()

    [rows, columns] = np.shape(R)
    if(rows!=columns):
        print("ERROR: Measurement noise matrix is not square.")
        sys.exit()

    [rows, columns] = np.shape(Pest0)
    if(rows!=columns):
        print("ERROR: Covariance matrix is not square.")
        sys.exit()

    if(len(xest0)!=rows):
        print("ERROR: State vector and covariance matrix are not compatible dimensions.")
        sys.exit()

    if not(np.allclose(Pest0, Pest0.T, rtol=1e-05, atol=1e-05)):
        print("ERROR: Covariance matrixx is not symmetric.")
        sys.exit()

    if not(callable(f)):
        print("ERROR: Input f must be callable.")
        sys.exit()

    if len(z)==1:
        [rows, columns] = np.shape(R)
        z = [z, [list(np.full((rows),np.nan))]]
    else:
        if h==[]:
            print("ERROR: Included measurements but missing measurement model.")
            sys.exit()

     # Estimating Dynamics Jacobian if not included based on time-marching scheme
    if(time == 'CT'):
        if(method == 'EE'):
            def df(f, x, dt, const):
                if const == []:
                    x1 = x + dt*f(x)
                else:
                    x1 = x + dt*f(x,const)
                return x1, dt
            if (F!=[]):
                def dF(f, df, F, x, dt, const):
                    if const == []:
                        return np.identity(len(x)) + dt*np.array(F(x))
                    else:
                        return np.identity(len(x)) + dt*np.array(F(x,const))
            else:
                def dF(f, df, F, x, dt, const):
                    dx = 1e-8
                    n = len(x)
                    df0, dx0 = df(f, x, dt, const)
                    jac = np.zeros((n, n))
                    for j in range(n):  
                        Dxj = (abs(x[j])*dx if x[j] != 0 else dx)
                        x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
                        dfi, dti = df(f,x_plus, dt, const)
                        jac[:, j] = (dfi - df0)/Dxj
                    return jac
        elif(method == 'RK4'):
            def df(f, x, dt, const):
                if const == []:
                    f1 = f(x)
                    f2 = f(x+(dt/2)*f1)
                    f3 = f(x+(dt/2)*f2)
                    f4 = f(x+dt*f3)
                else:
                    f1 = f(x,const)
                    f2 = f(x+(dt/2)*f1,const)
                    f3 = f(x+(dt/2)*f2,const)
                    f4 = f(x+dt*f3,const)
                x1 = x + dt*((f1/6)+(f2/3)+(f3/3)+(f4/6))
                return x1, dt
            def dF(f, df, F, x, dt, const):
                dx = 1e-8
                n = len(x)
                df0, dx0 = df(f, x, dt, const)
                jac = np.zeros((n, n))
                for j in range(n):  
                    Dxj = (abs(x[j])*dx if x[j] != 0 else dx)
                    x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
                    dfi, dti = df(f,x_plus, dt, const)
                    jac[:, j] = (dfi - df0)/Dxj
                return jac
        elif(method == 'RK45'):
            print("ERROR: RK45 is not working yet.")
            sys.exit()
        else:
            print("ERROR: Invalid time-marching scheme.")
            sys.exit()
    elif(time == 'DT'):
        def df(f, x, dt, const):
            if const == []:
                x1 = f(x,dt)
            else:
                x1 = f(x,dt,const)
            return x1, dt
        if (F!=[]):
            def dF(f, df, F, x, dt, const):
                if const == []:
                    return F(x,dt)
                else:
                    return F(x,dt,const)
        else:
            def dF(f, df, F, x, dt, const):
                dx = 1e-8
                n = len(x)
                df0, dx0 = df(f, x, dt, const)
                jac = np.zeros((n, n))
                for j in range(n):  
                    Dxj = (abs(x[j])*dx if x[j] != 0 else dx)
                    x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
                    dfi, dti = df(f,x_plus, dt, const)
                    jac[:, j] = (dfi - df0)/Dxj
                return jac

    # Estimating Measurement Jacobian if not included
    if(H == [])and(h!=[]):
        def dH(h,H,x):
            dx=1e-8
            n = len(x)
            func = h(x)
            jac = np.zeros((n, n))
            for j in range(n):  
                Dxj = (abs(x[j])*dx if x[j] != 0 else dx)
                x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
                jac[:, j] = (h(x_plus) - func)/Dxj
            return jac
    else:
        def dH(h,H,x):
            return H(x)

    tz = z[0]
    xest = np.array([xest0])
    Pest = np.array([Pest0])
    t0 = 0
    tspan = [t0]
    for i, tk1 in enumerate(tz):
        zk1 = z[1][i] 
        t0 = tspan[-1]

        # Prediction
        dtx = dt
        while t0 < tk1:
            dtx = min(dtx, tk1-t0)
            t0 += dtx
            tspan.append(t0)
            x1, dtx = df(f,xest[-1],dtx,const)
            xest = np.concatenate((xest,[x1]),axis=0)                                                                                                                      # Estimated state      
            Pest = np.concatenate((Pest,np.array([np.linalg.multi_dot([dF(f,df,F,xest[-1],dtx,const),Pest[-1],np.transpose(dF(f,df,F,xest[-1],dtx,const))]) + Q])),axis=0) # Estimated covariance                                                                                                           

        # Correction
        if not(np.isnan(zk1[0])):
            zpred = h(xest[-1])                                                                                                                             # Predicted Measurement
            K     = np.matmul(Pest[-1],np.transpose(dH(h,H,xest[-1])),np.linalg.inv(np.matmul(dH(h,H,xest[-1]),Pest[-1],np.transpose(dH(h,H,xest[-1])))+R)) # Kalman Gain
            xest[-1] = xest[-1] + np.matmul(K,zk1-zpred)                                                                                                    # Estimated state
            Pest[-1] = np.matmul((np.identity(len(xest0))-np.matmul(K,dH(h,H,xest[-1]))),Pest[-1])                                                          # Estimated covariance

    return xest, Pest, tspan
###########################################################################################################################################
def plot_gaussian_ellipsoid(m, C, sd=1, p=[], ax=[]):

    # Plot Gaussian ellipsoids representative of a 2D/3D covariance matrix

    # Parameters:
    #     m:   Mean
    #     C:   Covariance
    #     sd:  Number of standard deviations 
    #     ax:  Figure axis 
    #     p:   Plot parameters (optional)
    #         -type:    'fill' for filled in ellipse, 'line' for line plot ellipse
    #         -display: print to legend (1) or not (else)
    #         -name:    If display==1, name added to legend
    #         -color:   If 'line' type, color of line
    #         -lw:      If 'line' type, linewidth of ellipse
    #         -ls:      If 'line' type, line style of ellipse
    #         - means:  Plot means
    # Outputs:
    #     Plots Gaussian covariance matrix as an ellipse

    if(p == []):
        class ellipse: # Plotting Gaussians Ellipsoids
            def __init__(self, display, name, color, lw, ls, means):
                self.display = display
                self.name = name
                self.color = color
                self.lw = lw
                self.ls = ls
                self.means = means
        if(len(m))==2:
            p = ellipse(0, [], (0, 0, 1), 1, '-', 0)
        else:
            p = ellipse(0, [], (0, 0, 1, 0.2), [], [], 0)

    if(len(m)==2):
        npts = 50
        tt = np.linspace(0,2*np.pi, npts)
        x = np.cos(tt)
        y = np.sin(tt)
        ap = np.vstack([x,y])
        d, v = np.linalg.eig(C)
        d = np.diag(d)
        v = sd * np.sqrt(v)
        bp = np.linalg.multi_dot([v,d,ap])+np.transpose(np.tile(m,(np.size(ap[1]),1)))
        if(p.display):
            ax.plot(bp[0,:], bp[1,:], color = p.color, lw = p.lw, ls = p.ls, label = p.name)
        else:
            ax.plot(bp[0,:], bp[1,:], color = p.color)
        if(p.means):
            ax.scatter(m[0], m[1], s = 50,  marker = 'o', color = p.color)
    elif(len(m)==3):
        npts = 20
        theta = np.linspace(-npts,npts,(npts+1))/npts*np.pi
        phi = np.linspace(-npts,npts,(npts+1))/npts*np.pi/2
        cosphi = np.cos(phi)
        cosphi[0] = 0
        cosphi[npts] = 0
        sintheta = np.sin(theta)
        sintheta[0] = 0
        sintheta[npts] = 0
        x = np.outer(cosphi,np.cos(theta))
        y = np.outer(cosphi, sintheta)
        z = np.outer(np.sin(phi),np.full((1,npts+1),1))
        ap = np.vstack([np.ndarray.flatten(np.transpose(x)),np.ndarray.flatten(np.transpose(y)),np.ndarray.flatten(np.transpose(z))])
        d, v = np.linalg.eig(C)
        d = np.diag(d)
        d = sd * np.sqrt(d)
        bp = np.linalg.multi_dot([v,d,ap])+np.transpose(np.tile(m,(np.size(ap[1]),1)))
        xp = np.transpose(np.reshape(bp[0,:], np.shape(x)))
        yp = np.transpose(np.reshape(bp[1,:], np.shape(y)))
        zp = np.transpose(np.reshape(bp[2,:], np.shape(z)))
        ax.plot_surface(xp, yp, zp, color=p.color)
        if(p.means):
            ax.scatter(m[0], m[1], m[2], s = 50, marker = 'o', color = p.color)
    else:
        print("PLOT_GAUSSIAN_ELLIPSOIDS requires 2D or 3D means/covariances.")
        sys.exit()

    return
###########################################################################################################################################
def get_measurements(x0, dt, t, f, h, Q, R, method = 'RK4', const = []):

    # Given a dynamics model, initial state, and time parameters, return the
    # true states and measurements over time.
     
    # Inputs: 
    #    x0:      Initial true state (required)
    #    dt:      Time step (required)
    #    t:       Time parameter (required)
    #        - Time parameter can be in the following two fortmats:
    #            * t = [T] where T is period
    #            * t = [t1, t2, ..., T] where ti are the epochs
    #    f:       Dynamics model (required, continuous-time, function handle)
    #    h:       Measurement model (required, function handle)
    #    Q:       Process noise matrix (required)
    #    R:       Measurement noise matrix (required)
    #    method:  Time-marching method (optional)
    #        - 'EE':    Explicit Euler
    #        - 'RK4':   Runge-Kutta 4 (default)
    #        - 'RK45':  Adaptive Runge-Kutta 4/5
    #    const:   Miscellaneous constants (optional)

    #  Outputs:
    #    x:       True states
    #    z:       Measurements
    #    tx: True time span
    #    tz: Measurement time span

    if(method == 'EE'):
        def df(f, x, dt, const):
            x1 = x + dt*f(x,const)
            return x1, dt
    elif(method == 'RK4'):
        def df(f, x, dt, const):
            f1 = f(x,const)
            f2 = f(x+(dt/2)*f1,const)
            f3 = f(x+(dt/2)*f2,const)
            f4 = f(x+dt*f3,const)
            x1 = x + dt*((f1/6)+(f2/3)+(f3/3)+(f4/6))
            return x1, dt
    elif(method == 'RK45'):
        print("ERROR: RK45 is not working yet.")
        sys.exit()
    else:
        print("ERROR: Invalid time-marching scheme.")
        sys.exit()

    dtz = dt
    if len(t)==1:
        T = t[0]
        t = []
        t0 = 0
        while (t0 < T):
            t0 += dtz
            t.append(t0)
            dtz = min(dtz, T-t0)

    tx = [0]
    x = np.array([x0])
    tz = t
    t0 = 0
    cx = 0
    for i, tk1 in enumerate(tz):

        dtx = dt
        # True State
        while t0 < tk1:
            dtx = min(dtx, tk1-t0)
            t0 += dtx
            tx.append(t0)
            w  = np.linalg.multi_dot([np.sqrt(Q),np.random.randn(len(x0))])
            x1, dtx = df(f,x[cx,:],dtx,const)
            x = np.concatenate((x,[x1+w]),axis=0) 
            cx += 1

        # Measurement 
        d = R.shape[0]
        v = np.linalg.multi_dot([np.sqrt(R),np.random.randn(d)])
        if i==0:
            z = [h(x[cx,:]) + v]
        else:
            z = np.concatenate((z,[h(x[cx,:]) + v]),axis=0)

    return np.array(x), np.array(z), np.array(tx), np.array(tz)
###########################################################################################################################################
