import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import math

# Define constants and variables
Define_pipe_OD = np.array([1.5, 1.25, 1, 0.75, 0.5, 0.25, 0.125])
Define_pipe_ID = np.array([1.334, 1.084, 0.902, 0.694, 0.444, 0.194, 0.069])  # inches
Define_pipe_ID_m = Define_pipe_ID * 0.0254  # Convert ID from inches to meters

temperature_range = np.linspace(273.15, 305, 100)  # Kelvin
Substance = "N2O"
Mass_Flow = np.array([0.25, 0.5, 0.75, 1])  # kg/s
Dynamic_Viscosity = 1.47e-5  # Pa.s
Roughness = 0.015e-3  # meters
id_index = 4  # Corresponds to ID = 0.444 inches

# Create figures for plots
fig1, ax1 = plt.subplots(figsize=(8, 6))
fig2, ax2 = plt.subplots(figsize=(8, 6))
fig3, ax3 = plt.subplots(figsize=(8, 6))

# Find the index for 293.15 K in temperature_range
target_temp = 293.15
target_index = np.argmin(np.abs(temperature_range - target_temp))

# Perform calculations and plotting for each mass flow rate
for mass_flow in Mass_Flow:
    delta_p_l = np.zeros_like(temperature_range)
    flow_velocity = np.zeros_like(temperature_range)
    reynolds_numbers = np.zeros_like(temperature_range)
    densities = np.zeros_like(temperature_range)

    for i, temp in enumerate(temperature_range):
        density = CP.PropsSI('D', 'T', temp, 'Q', 0, Substance)
        densities[i] = density
        flow_velocity[i] = mass_flow / (density * np.pi * (Define_pipe_ID_m[id_index] ** 2) / 4)
        reynolds_num = (density * flow_velocity[i] * Define_pipe_ID_m[id_index]) / Dynamic_Viscosity
        reynolds_numbers[i] = reynolds_num
        relative_roughness = Roughness / Define_pipe_ID_m[id_index]

        if reynolds_num < 2300:
            friction_factor = 64 / reynolds_num
        else:
            theta1 = (-2.457 * np.log((7 / reynolds_num) ** 0.9 + (0.27 * relative_roughness))) ** 16
            theta2 = (37530 / reynolds_num) ** 16
            friction_factor = 8 * ((8 / reynolds_num) ** 12 + 1 / (theta1 + theta2) ** 1.5) ** (1 / 12)

        delta_p_l[i] = friction_factor * (density / 2) * (flow_velocity[i] ** 2) / Define_pipe_ID_m[id_index]

    # Plot results for current mass flow rate and label at 293.15 K
    ax1.plot(temperature_range, delta_p_l, label=f'Mass Flow: {mass_flow} kg/s')
    ax1.annotate(f'{delta_p_l[target_index]:.2f}', (target_temp, delta_p_l[target_index]), textcoords="offset points", xytext=(0,10), ha='center')

    ax2.plot(temperature_range, flow_velocity, label=f'Mass Flow: {mass_flow} kg/s')
    ax2.annotate(f'{flow_velocity[target_index]:.2f}', (target_temp, flow_velocity[target_index]), textcoords="offset points", xytext=(0,10), ha='center')

    ax3.plot(temperature_range, reynolds_numbers, label=f'Mass Flow: {mass_flow} kg/s')
    ax3.annotate(f'{reynolds_numbers[target_index]:.2f}', (target_temp, reynolds_numbers[target_index]), textcoords="offset points", xytext=(0,10), ha='center')

# Configure and display the plots
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('ΔP/L (Pa/m)')
ax1.set_title('Pressure Drop per Unit Length vs Temperature (.25" pipe)')
ax1.legend()

ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Flow Velocity (m/s)')
ax2.set_title('Flow Velocity vs Temperature (.25" pipe)')
ax2.legend()

ax3.set_xlabel('Temperature (K)')
ax3.set_ylabel('Reynolds Number')
ax3.set_title('Reynolds Number vs Temperature (.25" pipe)')
ax3.legend()

plt.show()
