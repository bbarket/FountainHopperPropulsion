from rocketcea.cea_obj import CEA_Obj,add_new_fuel,add_new_card,add_new_oxidizer
import math
import numpy as np
import matplotlib.pyplot as plt


# Define Fuel and Oxidizer parameters
#assume paraffin is only C25H52
card_str = """
fuel C25H52(S)   C 25 H 52   wt%=100
h,cal=-123000.0     t(k)=298.150
"""
N2O_str = """
oxid N2O(L) N 2 O 2 wt% = 100
h,cal = 19610.421 t(k) = 298.150
"""

# Add Fuel and Oxidizer to CEA
add_new_fuel( 'Paraffin', card_str )
add_new_oxidizer('NitrousOxide',N2O_str)

# Define CEA Object
ceaObject = CEA_Obj( oxName='NitrousOxide', fuelName='Paraffin')

#chamber pressure units are psia, eps = expansion ratio, frozen = 0 means we are in shifting equilibrium
# CEA Inputs
# Chamber pressure input is 600 psi (typical in range for hybrid rocket motors)
P0_CC_psi = 600
#Oxidizer to Fuel Mass Ratio - defined parameter from graphs for Isp and C* for Paraffin - N20, 7.5 is known good
OF = 7.5
# Nozzle Expansion Ratio AR = 4, range is typically from 3-25, but we are operating with smaller motor so more like 3-15
AR = 4

#parameters for get full cea output
#Pc=100.0, MR=1.0, eps=40.0, subar=None, PcOvPe=None, frozen=0, frozenAtThroat=0, short_output=0, show_transport=1, pc_units='psia', output='calories', show_mass_frac=False, fac_CR=None
s = ceaObject.get_full_cea_output( Pc=600, MR=7.5, eps=4, short_output=1,frozen=0)
print(s)

'''
c_star_metric = ceaObject.get_Cstar(Pc = P0_CC_psi, MR = OF) * .3048
c_star_imperial = ceaObject.get_Cstar(Pc = P0_CC_psi, MR = OF)
print("c* in m/s",c_star_metric)

gamma_combustionChamber = ceaObject.get_Chamber_MolWt_gamma(Pc = P0_CC_psi,MR = OF, eps = AR)[1]
print("Gamma in combustion chamber",gamma_combustionChamber)

pressure_exit_psia = 1/(1/P0_CC_psi* ceaObject.get_PcOvPe(Pc=P0_CC_psi,MR = OF, eps = AR, frozen = 0, frozenAtThroat = 0))
pressure_exit_pa = pressure_exit_psia*6894.76
print("Exit pressure (Pa)",pressure_exit_pa)

C_f = ceaObject.get_PambCf(Pc = P0_CC_psi,MR = OF, eps = AR)[1]
print("C_f (nozzle efficiency)",C_f)

'''

def getExitPressureATM(chamberPressure,OF,AR):
    pressure_exit_psia = 1/(1/simChamberPressure* ceaObject.get_PcOvPe(chamberPressure,OF, AR, frozen = 0, frozenAtThroat = 0))
    pressure_exit_pa = pressure_exit_psia*6894.76
    pressure_exit_atm = pressure_exit_pa/101325
    return pressure_exit_atm









pi = math.pi
'''

#----------------------------------------------------------------------------
# Outputs and Key Parameters
# now from CEA results we will calculate key design parameters using formulas from Prof. Cantwell
pi = math.pi
#g_imperial is 32.174 ft/s^2
g_imperial = 32.174
#Thrust = 600 lbf
Thrust = 600
#burn time is in seconds
burn_time = 15
#Max OD Grain (inches) = 4.875 assuming liner thickness of 1/16"
OD_grain = 4.875
OD_grain_metric = OD_grain * .0254
# Grain max radius
OR_grain_metric = OD_grain_metric/2
#simulation constants (from Cantwell)
a0 = 1.55*10**-4
n = .5
#Initial Oxidizer Tank Pressure (Pa)  (1000 psi)
P0_ox = 6.895*10**6
#density of paraffin fuel (kg/m^3)
rho_fuel = 924.5
# Initial Thrust Target (N) (600 lbf)
F0 = 2668.93
# Initial Chamber Pressure (600 psi)
P0_CC = 4.137*10**6
#Nozzle Area Ratio (dimensionless)
#previously defined
#Efficiency of C* (typically between 75 and 90% for small hybrids)
n_c = .85
#Discharge Coefficient (assumed to be one for first pass calculations)
C_d = 1
#Nozzle Efficiency (takes into account loss due to conical nozle (1.7%), loss due to boundar layer (1%), loss due to real gas props (.4%), loss due to small particles (2%))
n_nozzle_divergence = .983
n_nozzle_friction = .99
n_nozzle_realGasProps = .996
n_nozzle_smallParticles = .98
n_n = n_nozzle_divergence*n_nozzle_friction*n_nozzle_realGasProps*n_nozzle_smallParticles
#Nozzle Half Angle (should be 15d for conical nozzle, reccomended from Hybrid Handbook)
nozzle_half_angle = 15
'''
#From CEA/Inputs
#Nozzle sizing for expansion section (need to identify sizing params for throat length and subsonic section)(units will be inches)
def getThroatDistances(Thrust,chamberPressure,nozzleEfficiency,C_f,dischargeCoefficient,nozzleHalfAngle,AreaRatio):
    # Nozzle Throat Area -  Equation 7.14 Hybrid Handbook A_throat = F_thrust/ (C_f * C_d * n_nozzle * P_chamber)
    # units = lbf / (lbf / in^2 * unitless * unitless) = sq inches
    areaThroatInchSquare = Thrust / (chamberPressure * nozzleEfficiency* C_f * dischargeCoefficient)
    radiusThroatInch = (areaThroatInchSquare / (pi))**.5
    radiusExitInch = AreaRatio**1/2 * radiusThroatInch
    lengthThroatToExit = (radiusExitInch-radiusThroatInch)/(math.tan(math.radians(nozzleHalfAngle)))
    return(radiusThroatInch,radiusExitInch,lengthThroatToExit)

# Mass flow rate of propellant (eqn 7.15 from hybrid handbook) (units will be lb/sec and kg/sec returned
def getMassFlowPropellant(chamberPressure,dischargeCoefficient,throatArea,cStarImperial,cStarEfficiency):
    massFlowPropellantImperial = gImperial*chamberPressure * dischargeCoefficient * throatArea/ (cStarImperial*cStarEfficiency)
    massFlowPropellantMetric = massFlowPropellantImperial*0.453592
    return massFlowPropellantImperial,massFlowPropellantMetric

# With mass flow rate of propellant and OF ratio we can get mass flow of N2O and Paraffin
def getMassFlowFuel(massFlowPropellant, OF):
    massFlowFuelImperial = massFlowPropellant/(OF+1)
    massFlowFuelMetric = massFlowFuelImperial * 0.453592
    return massFlowFuelImperial,massFlowFuelMetric
def getMassFlowOxidizer(massFlowPropellant,OF):
    massFlowOxImperial = OF * massFlowPropellant / (OF+1)
    massFlowOxMetric = massFlowOxImperial * 0.453592
    return massFlowOxImperial,massFlowOxMetric
def getTotalOx(massFlowOxImperial,burnTime):
    massTotalImperialOx = massFlowOxImperial*burnTime
    massTotalMetricOx = massTotalImperialOx * 0.453592
    return massTotalImperialOx,massTotalMetricOx
def getTotalFuel(massFlowFuelImperial,burnTime):
    massTotalImperialFuel = massFlowFuelImperial*burnTime
    massTotalMetricFuel = massTotalImperialFuel * 0.453592
    return massTotalImperialFuel,massTotalMetricFuel
# To find initial grain diameter we can use eqn 7.21 plugging in max time of burn and max diameter to solve for d_initial, DENSITY IS IN INCHES^3
def getGrainSizing(outerDiameterGrainImperial,a0Coefficient,nCoefficient,burnTime,massFlowOxImperial,massFuelImperial,densityFuelImperial):
    massFlowOxMetric = massFlowOxImperial * 0.453592
    outerDiameterGrainMetric = outerDiameterGrainImperial * .0254
    diameterInitialMetric = (outerDiameterGrainMetric ** (2*nCoefficient+1) - (2*nCoefficient+1)*2**(2*nCoefficient+1)*a0Coefficient*massFlowOxMetric**nCoefficient*burnTime / pi**nCoefficient)**(1/(2*nCoefficient+1))
    
    diameterInitialImperial = diameterInitialMetric/.0254
    volumeFuelMetric = massFuelImperial*0.453592 / 924.5
    volumeFuelImperial = massFuelImperial/densityFuelImperial


    #print("outerDiameterGrainImperial:", outerDiameterGrainImperial)
    #print("diameterInitialImperial:", diameterInitialImperial)
    #print("volumeFuelImperial:", volumeFuelImperial)
    #print("Calculated Area:", (pi / 4 * (outerDiameterGrainImperial ** 2 - diameterInitialImperial ** 2)))

    lengthGrainImperial = volumeFuelImperial / (pi/4*(outerDiameterGrainImperial**2 - diameterInitialImperial**2))
    
    return diameterInitialImperial,volumeFuelImperial,lengthGrainImperial


# From CEA we need to get C_f,C*
#Simulation Conditions: Chamber Pressure = 600 psia, Thrust = 600lbf, AR = 1-40: 100 total, OF = 7.5
#To do sizing we only need C_f and C*
# We want to see throat/exit sizes, mass flow rates, and total masses of fuel and oxidizer as f(AR)
simThrust = 600
simBurnTime = 15
simDischargeCoefficient = 1
#simNozzleEfficiency = .94989
simNozzleEfficiency = 1
simCStarEfficiency = .85
simHalfAngle = 15
simChamberPressure = 600
simOuterDiameterGrain = 4.875
a0=1.55*10**-4
n=.5
#density is lb/in^3
densityFuelImperial = 0.0334
simOF = 7.5
simAR = np.linspace(2,15,100)
cF = np.linspace(1,40,100)
cStar = np.linspace(1,40,100)
gImperial = 32.174
# Initialize lists to store results
throat_lengths = []
exit_lengths = []
nozzle_lengths = []
mass_flow_oxidizers_metric = []
total_mass_oxidizers_metric = []
total_mass_fuels_metric = []
grain_length = []
grain_port_diameter = []
exit_pressure = []

count = 0
for i in simAR:
    cF[count] = ceaObject.get_PambCf(Pc = simChamberPressure,MR = simOF, eps = i)[1]
    cStar[count] = ceaObject.get_Cstar(Pc = simChamberPressure, MR = simOF)
    exit_pressure.append(getExitPressureATM(simChamberPressure,simOF,i))
    
    throatRad,exitRad,nozzleLength = getThroatDistances(simThrust,simChamberPressure,simNozzleEfficiency,cF[count],simDischargeCoefficient,simHalfAngle,i)
    massFlowPropImperial,massFlowPropMetric = getMassFlowPropellant(simChamberPressure,simDischargeCoefficient,pi*throatRad**2,cStar[count],simCStarEfficiency)
    massFlowOxImperial,massFlowOxMetric = getMassFlowOxidizer(massFlowPropImperial,simOF)
    massFlowFuelImperial,massFlowFuelMetric = getMassFlowFuel(massFlowPropImperial,simOF)
    totalMassOxImperial,totalMassOxMetric = getTotalOx(massFlowOxImperial,simBurnTime)
    totalMassFuelImperial,totalMassFuelMetric = getTotalFuel(massFlowFuelImperial,simBurnTime)
    portDiameterImperial,volumeFuelImperial,lengthGrainImperial = getGrainSizing(simOuterDiameterGrain,a0,n,simBurnTime,massFlowOxImperial,totalMassFuelImperial,densityFuelImperial)

    throat_lengths.append(throatRad)  # or throatRad if you want the radius
    exit_lengths.append(exitRad)
    nozzle_lengths.append(nozzleLength)
    mass_flow_oxidizers_metric.append(massFlowOxMetric)
    total_mass_oxidizers_metric.append(totalMassOxMetric)
    total_mass_fuels_metric.append(totalMassFuelMetric)
    grain_length.append(lengthGrainImperial)
    #print(lengthGrainImperial)
    grain_port_diameter.append(portDiameterImperial)
    

    count +=1
    
"""
# Example: Plotting nozzle lengths vs AR
plt.subplot(3, 2, 1)
plt.plot(simAR, nozzle_lengths, marker='o', label='Nozzle Length')
plt.xlabel('Area Ratio (AR)')
plt.ylabel('Nozzle Length (units)')
plt.title('Nozzle Length vs AR')
plt.grid()
plt.legend()
"""
# Example: Plotting total oxidizer mass vs AR
plt.plot(simAR, mass_flow_oxidizers_metric, linestyle='-', label='Total Oxidizer Mass Flow', color='orange')
plt.xlabel('Area Ratio (AR)')
plt.ylabel('Total Oxidizer Mass Flow (kg/s)')
plt.title('Total Oxidizer Mass Flow vs AR')
plt.grid()
plt.legend()

# Additional plots for other metrics can be added similarly
plt.show()

plt.plot(simAR, total_mass_oxidizers_metric, linestyle='-', label='Total Oxidizer Mass', color='orange')
plt.xlabel('Area Ratio (AR)')
plt.ylabel('Total Oxidizer Mass (kg)')
plt.title('Total Oxidizer Mass vs AR')
plt.grid()
plt.legend()
plt.show()

plt.plot(simAR, exit_lengths, linestyle='-', label='Exit Radius (inch)', color='orange')
plt.xlabel('Area Ratio (AR)')
plt.ylabel('Exit Radius (inch)')
plt.title('Exit Radius vs AR')
plt.grid()
plt.legend()
plt.show()

plt.plot(simAR, grain_port_diameter, linestyle='-', label='Grain Port Diameter (inch)', color='orange')
plt.xlabel('Area Ratio (AR)')
plt.ylabel('Grain Port Diameter (inch)')
plt.title('Grain Port Diameter vs AR')
plt.grid()
plt.legend()
plt.show()

plt.plot(simAR, grain_length, linestyle='-', label='Grain Length (inch)', color='orange')
plt.xlabel('Area Ratio (AR)')
plt.ylabel('Grain Length (inch)')
plt.title('Grain Length vs AR')
plt.grid()
plt.legend()
plt.show()
plt.plot(simAR, cF, linestyle='-', label='C_f (nozzle coefficient)', color='orange')
plt.xlabel('Area Ratio (AR)')
plt.ylabel('C_f')
plt.title('C_f vs AR')
plt.grid()
plt.legend()
plt.show()

plt.plot(simAR, exit_pressure, linestyle='-', label='Exit Pressure (atm)', color='orange')
plt.xlabel('Exit Pressure (atm)')
plt.ylabel('Exit Pressure')
plt.title('Exit Pressure vs AR')
plt.grid()
plt.legend()
plt.show()
















