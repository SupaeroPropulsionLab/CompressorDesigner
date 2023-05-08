import math

# backswept impeller

# Known design inlet conditions

absolute_inlet_stagnation_temp = 288 # K TT1
absolute_inlet_stagnation_pres = 1.01e5 # Pa
c1 = 150 # m/s
M1 = 0.5 # [1]
up_to_tip_radius_ratio = 0.5 # [1]
Meye = 1 # [1]
compressor_adiabatic_efficiency = 0.72 # [1]
mass_flow_rate = 5 # kg/s
cp = 1004 
gamma = 1.4

# Known exit conditions

cr2 = 150 # m/s
cz1 = 150 # m/s
u2 = 520 # m/s
beta2_prime = 20 # deg
epsilon = 10 # deg
beta2 = 30 # deg

def static_temp_inlet(Tt1,c2,cp):
    return Tt1 - c2**2/(2*cp)

def a1(T1):
    return math.sqrt((gamma - 1)*cp*T1)

def M1(c1,a1):
    return c1/a1

def C_theta2(u2,c1,beta2):
    return u2 - c1*math.tan(beta2*math.pi/180)

# absolute velocity at exit
def c2(c_theta2,cr2):
    return math.sqrt(c_theta2**2 + cr2**2)

def a2(a1,MT2,M1):
    return a1*math.sqrt(1 + ((gamma - 1)/2 * MT2**2)/(1 + (gamma - 1)/2 * M1**2))

def Tt2(Tt1, u2, c_theta2):
    return Tt1*(1 + u2**2/(cp*Tt1)*(c_theta2/u2))

def T2(Tt2,c2):
    return Tt2 - c2**2/(2*cp)
    
T1 = static_temp_inlet(absolute_inlet_stagnation_temp,cr2,cp)
a_1 = a1(T1)
M1 = M1(c1,a_1)
c_theta2 = C_theta2(u2,c1,beta2)
c_2 = c2(c_theta2,cr2)
a_2 = a2(a_1, u2/a_1, M1)
T_t2 = Tt2(absolute_inlet_stagnation_temp, u2, c_theta2)
T_2 = T2(T_t2,c_2)
a_2_exact = a1(T_2)
# impeller exit mach number
M_theta2 = c_theta2/a_2
# M_r2 = cr2/a_2
M2 = c_2/a_2
M_r2 = (cr2/math.cos(beta2*math.pi/180))/a_2
tau_c = T_t2/absolute_inlet_stagnation_temp 
pi_c = (1 + compressor_adiabatic_efficiency*(tau_c - 1))**(gamma/(gamma - 1))
p_t2 = absolute_inlet_stagnation_pres*pi_c
# compressor polytropic efficiency
e_c = ((gamma - 1)/gamma)*math.log(pi_c)/math.log(1 + (pi_c**((gamma - 1)/gamma) - 1)/compressor_adiabatic_efficiency)
# sizing of inlet flow
p1 = absolute_inlet_stagnation_pres/((1 + (gamma - 1)/2*M1**2)**((gamma)/(gamma - 1)))
print(p1)
rho_1 = gamma*p1/(a_1**2)
A1 = mass_flow_rate/(rho_1*cz1)