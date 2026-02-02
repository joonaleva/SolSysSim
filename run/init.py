import numpy as np
import csv
import os

## CONSTANTS ##
G_SI = 6.67430e-11  #m^3 kg^-1 s^-2
DAY = 86400  #s
YEAR = 365.25*DAY  #s, Julian year
AU = 1.495978707e11  #m
PC = 206.265*AU  #m
M_SOL = 1.988416e30  #kg
M_EARTH = 5.972168e24 #kg
G_SOL = G_SI / (PC**3 * M_SOL**(-1) * YEAR**(-2)) #pc^3 M_sol^-1 year^-2
G_AUD = G_SI / (AU**3*DAY**(-2))  #au^3 kg^-1 day^-2

## SIMULATION CONFIG ## (USE INTEGERS)
Units = 2                           # 1:m|s|kg  2:au|day|kg  3:pc|yr|Msol
Center = 0                          # Index of center object,  0: center of mass
Duration = 5000                    # Simulation duration in chosen units
Steps = int(Duration*100)           # No. of steps, currently 100 steps per day
Print_step = int(Steps/10)          # Interval of console print outs, currently 10 times overall
Save_step = int(Steps/Duration/10)  # Interval of position saves into output data file, currently ten times a day

## SETTINGS ##
CurrTime = 2454194          # Current time in JD, 19.12.2025 = 2461029
IncludeDwarfs = True        # Include dwarf planets
IncludeMoons = True         # Include moons

Simulation = False           # Run the simulation directly
Plot = False                 # Do plot directly

## TARGET LIST ##
                            # a=semi-major axis, e=eccentricity, i=inclination, Ω=longitude of ascending node, 
                            # ω=argument of periapsis, t_peri=time of periapsis, m=mass         (all angles in degrees)
planets = [
    # name, a, e, i, Ω, ω, t_peri, m                #Orbits are around the first object in the list
    ("Sun", 0, 0, 0, 0, 0, CurrTime, M_SOL),
    ("Mercury", 0.387098, 0.205630, 7.005, 48.331, 29.124, 2460740, 3.3011e23),   
    ("Venus", 0.723332, 0.006772, 3.39458, 76.680, 54.884, 2460717, 4.8675e24),
    ("Earth", 1, 0.0167086, 0, -11.26064, 114.20783, 2460680, M_EARTH),
    ("Mars", 1.52368055, 0.0934, 1.850, 49.57854, 286.5, 2460439, 5.972168e24),
    ("Jupiter", 5.2038, 0.0489, 1.303, 100.464, 273.867, 2459966, 1.8982e27),
    ("Saturn", 9.5826, 0.0565, 2.485, 113.665, 339.392, 2463566, 5.6834e26),
    ("Uranus", 19.19126, 0.04717, 0.773, 74.006, 96.998857, 2470037, 8.6810e25),
    ("Neptune", 30.07, 0.008678, 1.770, 131.783, 273.187, 2467132, 1.02409e26),
    #("Vulcan", 0.14, 0.04, 12.167, 45, 50, CurrTime, 8e22),                     #Objects can be commented out
    #("Planet Nine", 290, 0.29, 6.8, 90, 180, CurrTime, 4.4*M_EARTH),
]

dwarfs = [
    # name, a, e, i, Ω, ω, t_peri, m                    #Orbits are around the first planet
    ("Ceres", 2.77, 0.0785, 10.6, 80.3, 73.6, 2459921, 9.3839e20),
    ("Pluto", 39.482, 0.2488, 17.16, 110.299, 113.834, 2447775, 1.3025e22),
    ("Eris", 67.99, 0.4370, 43.86, 36.02, 150.73, 2545752, 1.638e22),
    #("Sedna", 506, 0.8496, 11.9307, 144.248, 311.352, 2479503, 1e21),
    ("Oumuamua", -1.2723, 1.20113, 122.74, 24.597, 241.811, 2458006, 1e10),
]

moons = [
    # name, a, e, i, m, index of parent planet (note python indexes, Sun is 0, Mercury is 1 etc)
    ("Moon", 1/389, 0.0549, 5.145, 7.346e22, 3),
    #("Phobos", 9.376e6/AU, 0.0151, 26.04, 1.060e16, 4),                # Objects can be commented out
    #("Deimos", 23.4632e6/AU, 0.00033, 27.58, 1.51e15, 4),
    ("Io", 0.4127e9/AU, 0.004031, 2.213, 8.931938e22, 5),
    ("Europa", 0.6709e9/AU, 0.009, 1.791, 4.79984e22, 5),
    ("Ganymede", 1.0704e9/AU, 0.0013, 2.214, 1.4819e23, 5),
    ("Callisto", 1.8827e9/AU, 0.0074, 2.017, 1.075938e23, 5),
    ("Enceladus", 0.238037e9/AU, 0.0047, 0, 1.080318e20, 6),           
    ("Titan", 1.221870e9/AU, 0.0288, 0, 1.34518e23, 6),
    ("Triton", 0.354759e9/AU, 0.000016, 129.812, 2.1389e22, 8),
]

dwarfmoons = [
    # name, a, e, i, m, parent planet (note python indexes)
    ("Charon", 1.9595764e7/AU, 0.000161, 112.783, 1.5897e21, len(planets)+1),
    ("Dysnomia", 0.037273e9/AU, 0.0062, 61.59, 8.2e19, len(planets)+2),
]



############################# ACTUAL CODE BEGINS, DO NOT EDIT FURTHER FROM HERE ####################################

# count target object amount
if IncludeDwarfs: planets = planets+dwarfs
N_obj = len(planets)
if IncludeMoons: 
    if IncludeDwarfs: moons = moons+dwarfmoons
    N_obj += len(moons)

## POSITION AND VELOCITY AT PERIHELION CALCULATOR ##
def rv_at_perihelion(a, e, i_deg, Omega_deg, omega_deg, mu):
    #semi major axis, eccentricity, inclination, Lgtude of asc node, Arg of peri, Grav parameter

    #failsafe for a=0
    if a==0:
        r_inertial = np.array([0,0,0])
        v_inertial = np.array([0,0,0])
        return r_inertial, v_inertial

    # convert to radians
    i = np.deg2rad(i_deg)
    Omega = np.deg2rad(Omega_deg)
    omega = np.deg2rad(omega_deg)

    # perihelion distance
    r_p = a * (1.0 - e)

    # specific angular momentum
    h = np.sqrt(mu * a * (1.0 - e**2))

    # perifocal vectors (f=0)
    r_pf = np.array([r_p, 0.0, 0.0])
    v_pf = np.array([0.0, mu/h * (1.0 + e), 0.0])  # or use v_p = sqrt(mu*(1+e)/(a*(1-e)))

    # rotation matrices
    def Rz(theta):
        c,s = np.cos(theta), np.sin(theta)
        return np.array([[c,-s,0],[s,c,0],[0,0,1]])
    def Rx(phi):
        c,s = np.cos(phi), np.sin(phi)
        return np.array([[1,0,0],[0,c,-s],[0,s,c]])

    # perifocal to equatorial conversion
    R = Rz(Omega) @ Rx(i) @ Rz(omega)

    r_inertial = R @ r_pf
    v_inertial = R @ v_pf

    return r_inertial, v_inertial


## POSTITION AND VELOCITY AT ANY TIME ##
def rv_at_time(a, e, i, Omega, omega, t_peri, mu):

    if a==0:                            #failsafe for a=0
        r_vec = np.array([0,0,0])
        v_vec = np.array([0,0,0])
        return r_vec, v_vec
    
    # convert to radians
    i = np.deg2rad(i)
    Omega = np.deg2rad(Omega)
    omega = np.deg2rad(omega)

    # rotation matrices
    def Rz(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c,-s,0],[s,c,0],[0,0,1]])

    def Rx(phi):
        c, s = np.cos(phi), np.sin(phi)
        return np.array([[1,0,0],[0,c,-s],[0,s,c]])

    # mean motion (works for all conics except parabola)
    if e != 1.0:
        n = np.sqrt(mu /(abs(a)**3))     # with kepler 3 and n=2pi/P

    # time since periapsis
    dt = (CurrTime - t_peri)

    ## Solve true anomaly and heliocentric distance ##  (wikipedia: Kepler's laws of planetary motion)
    if e < 1:  
        # elliptic: E - e sin E = M
        M = n * dt                      #mean anomaly
        E = kepler_E(M, e)              #eccentric anomaly
        r = a * (1 - e*np.cos(E))       #heliocentric distance
        f = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))    #true anomaly

    elif e > 1:
        # hyperbolic: e sinh H - H = M
        M = -n * dt             #added minus cause the direcion in simu was wrong for some reason 
        H = kepler_H(M, e)      
        r = a * (e*np.cosh(H) - 1)
        f = 2*np.arctan2(np.sqrt(e+1)*np.sinh(H/2), np.sqrt(e-1)*np.cosh(H/2))

    else:
        # parabolic: use Barker's equation
        D = (3*np.sqrt(mu/(2*a*a*a)) * dt)**(1/3)
        r = a*(1 + D*D)
        f = 2*np.arctan(D)


    ## Perifocal position/velocity  ##
    h = np.sqrt(mu * a * (1 - e*e))
    r_pf = np.array([r*np.cos(f), r*np.sin(f), 0])

    v_pf = np.array([-mu/h * np.sin(f), mu/h * (e + np.cos(f)), 0])


    ## Rotate into inertial frame  ##       (perifocal to equatorial conversion)
    R = Rz(Omega) @ Rx(i) @ Rz(omega)

    r_vec = R @ r_pf
    v_vec = R @ v_pf

    return r_vec, v_vec


# NUMERICAL KEPLER EQUATION SOLVERS
def kepler_E(M, e, tol=1e-12):
    #Solve elliptical Kepler's equation E - e sin E = M.   #JPL: Approximate positions of the planets
    E = M
    for _ in range(50):
        dE = M -(E - e*np.sin(E)) / (1 - e*np.cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def kepler_H(M, e, tol=1e-12):
    #Solve hyperbolic Kepler's equation e sinh H - H = M.
    # initial guess
    H = np.log(2*M/e + 1.8) if M > 0 else -np.log(-2*M/e + 1.8)
    for _ in range(50):
        f = e*np.sinh(H) - H - M
        fp = e*np.cosh(H) - 1
        H -= f/fp
        if abs(f) < tol:
            break
    return H


### WRITE DATA COLUMNS INTO FILE ###
names = []
rx, ry, rz = [], [], []
vx, vy, vz = [], [], []
mass = []

if Units==1: G_para = G_SI
if Units==2: G_para = G_AUD
if Units==3: G_para = G_SOL

### FOR PLANETS ###
for name, a, e, inc, Om, om, t_peri, m in planets:
    r, v = rv_at_time(a, e, inc, Om, om, t_peri, G_para*(planets[0][-1]+m))
    names.append(name)
    rx.append(r[0]); ry.append(r[1]); rz.append(r[2])
    vx.append(v[0]); vy.append(v[1]); vz.append(v[2])
    mass.append(m)

### FOR MOONS ###
if IncludeMoons:
    for name, a, e, inc, m, parent in moons:
        r, v = rv_at_perihelion(a, e, inc, 0, 0, G_para*(planets[parent][-1]+m))
        names.append(name)
        rx.append(rx[parent]+r[0]); ry.append(ry[parent]+r[1]); rz.append(rz[parent]+r[2])
        vx.append(vx[parent]+v[0]); vy.append(vy[parent]+v[1]); vz.append(vz[parent]+v[2])
        mass.append(m)

# Write config file
with open("config.txt", "w", newline="") as f:
    writer = csv.writer(f, delimiter='\t')      #'\t' is tab and ' ' is space
    writer.writerow([N_obj, Units, Center, Duration, Steps, Print_step, Save_step])

# Write target sequence for plotting
with open("sequence.txt", "w", newline="") as f:
    writer = csv.writer(f, delimiter='\n')      #'\t' is tab, ' ' is space, '\n' is newline
    writer.writerow(names)

# Write data file in columns
with open("input.dat", "w", newline="") as f:
    writer = csv.writer(f, delimiter='\t')      #'\t' is tab and ' ' is space

    writer.writerow(mass)
    writer.writerow(rx)
    writer.writerow(ry)
    writer.writerow(rz)
    writer.writerow(vx)
    writer.writerow(vy)
    writer.writerow(vz)

if Simulation: os.system("./simulation")
if Plot: os.system("python3 plot.py")