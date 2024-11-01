# Date Created: Oct 31, 2024
# Author: Brandon Barket
# Purpose: Hosts functions for hybrid sizing and uses input from rocketCEA and design parameters, improves organization over previous document rocketCEATester.py
from rocketcea.cea_obj import CEA_Obj,add_new_fuel,add_new_card,add_new_oxidizer
import math
import numpy as np
import matplotlib.pyplot as plt

# Constants
pi = math.pi
gImperial = 32.174


#-------------- CEA SETUP --------------

#Define Fuel and Oxidizer
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

# CEA FUNCTIONS
#returns Cf from chamber pressure, OF, and AR
def getCf(ChamberPressure,OF,AR):
    # requires amb pressure, chamber pressure, OF and AR
    return ceaObject.get_PambCf(14.7,ChamberPressure,OF,AR)[0]

#return C* in ft/s from chamber pressure and OF
def getCStarFt(ChamberPressure,OF):
    return ceaObject.get_Cstar(ChamberPressure,OF)


# this will print traditional CEA output
#s = ceaObject.get_full_cea_output(600,7.5,4)



# -------------- Sizing Function Definitions --------------

# Sizes the nozzle using C_f and area ratio of nozzle
def getThroatDistances(Thrust,chamberPressure,nozzleEfficiency,C_f,dischargeCoefficient,nozzleHalfAngle,AreaRatio):
    # Nozzle Throat Area -  Equation 7.14 Hybrid Handbook A_throat = F_thrust/ (C_f * C_d * n_nozzle * P_chamber)
    # units = lbf / (lbf / in^2 * unitless * unitless) = sq inches
    areaThroatInchSquare = Thrust / (chamberPressure * nozzleEfficiency* C_f * dischargeCoefficient)
    radiusThroatInch = (areaThroatInchSquare / (pi))**.5
    radiusExitInch = AreaRatio**(1/2) * radiusThroatInch
    lengthThroatToExit = (radiusExitInch-radiusThroatInch)/(math.tan(math.radians(nozzleHalfAngle)))
    return(radiusThroatInch,radiusExitInch,lengthThroatToExit)

# Mass flow rate of propellant (eqn 7.15 from hybrid handbook) (units will be lb/sec and kg/sec returned
def getMassFlowPropellant(chamberPressure,dischargeCoefficient,throatArea,cStarImperial,cStarEfficiency):
    massFlowPropellantImperial = gImperial*chamberPressure * dischargeCoefficient * throatArea/ (cStarImperial*cStarEfficiency)
    massFlowPropellantMetric = massFlowPropellantImperial*0.453592
    return massFlowPropellantImperial,massFlowPropellantMetric

# With mass flow rate of propellant and OF ratio we can get mass flow of N2O and Paraffin and then total masses from burn time
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
    #volumeFuelMetric = massFuelImperial*0.453592 / 924.5
    volumeFuelImperial = massFuelImperial/densityFuelImperial


    #print("outerDiameterGrainImperial:", outerDiameterGrainImperial)
    #print("diameterInitialImperial:", diameterInitialImperial)
    #print("volumeFuelImperial:", volumeFuelImperial)
    #print("Calculated Area:", (pi / 4 * (outerDiameterGrainImperial ** 2 - diameterInitialImperial ** 2)))

    lengthGrainImperial = volumeFuelImperial / (pi/4*(outerDiameterGrainImperial**2 - diameterInitialImperial**2))
    
    return diameterInitialImperial,volumeFuelImperial,lengthGrainImperial

# -------------- Example CEA + Sizing --------------
simThrust = 600
simBurnTime = 15
simDischargeCoefficient = 1
simNozzleEfficiency = .94989
simCStarEfficiency = .85
simHalfAngle = 15
simChamberPressure = 600
simOuterDiameterGrain = 4.875
a0=1.55*10**-4
n=.5
#density of paraffin is in lb/in^3
densityFuelImperial = 0.0334
simOF = 7.5
simAR = 6.8

# Get Cf and C* from CEA

simCf = getCf(simChamberPressure,simOF,simAR)
simCStar = getCStarFt(simChamberPressure,simOF)


simThroatRad, simExitRadius,simThroatToExitLength = getThroatDistances(simThrust,simChamberPressure,simNozzleEfficiency,simCf,simDischargeCoefficient,simHalfAngle,simAR)
simThroatArea = pi*simThroatRad**2
simMassFlowPropImperial, simMassFlowPropMetric = getMassFlowPropellant(simChamberPressure,simDischargeCoefficient,simThroatArea,simCStar,simCStarEfficiency)
simMassFlowFuelImperial,simMassFlowFuelMetric = getMassFlowFuel(simMassFlowPropImperial,simOF)
simMassFlowOxImperial,simMassFlowOxMetric = getMassFlowOxidizer(simMassFlowPropImperial,simOF)
simTotalFuelImperial,simTotalFuelMetric = getTotalFuel(simMassFlowFuelImperial,simBurnTime)
simTotalOxidizerImperial,simTotalOxidizerMetric = getTotalOx(simMassFlowOxImperial,simBurnTime)
simGrainPortDiameter, simVolumeFuelImperial,simLengthGrain = getGrainSizing(simOuterDiameterGrain,a0,n,simBurnTime,simMassFlowOxImperial,simTotalFuelImperial,densityFuelImperial)
print("--- CEA VALUES ---")
print("simCf",simCf)
print("simC*",simCStar)

print("--- Calculated Values and Targets ---")
print("Thrust (lbf)",simThrust)
print("Burn Time (s)",simBurnTime)
print("Expansion Ratio",simAR)
print("Exit Pressure",1/(1/simChamberPressure* ceaObject.get_PcOvPe(simChamberPressure,simOF, simAR))*6894.76/101325)
print("Chamber Pressure (psia)",simChamberPressure)
print("Throat Diameter (in)",2*simThroatRad)
print("Grain Port Diameter (in)",simGrainPortDiameter)
print("Grain Length (in)",simLengthGrain)
print("Mass of Fuel (lb)",simTotalFuelImperial)
print("Mass of Oxidizer (lb)",simTotalOxidizerImperial)
print("Ratio of Grain OD / Grain Port Diameter",simOuterDiameterGrain/simGrainPortDiameter)
