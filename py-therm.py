from dataclasses import dataclass
from typing import List
import math
import matplotlib.pyplot as plt


def calculate_vapor_pressure(temperature_celsius: float, relative_humidity: float) -> float:
    """
    Calculate the vapor pressure at a given temperature and relative humidity.

    Parameters:
    temperature_celsius (float): Temperature in degrees Celsius.
    relative_humidity (float): Relative humidity as a percentage (0-100).

    Returns:
    float: Vapor pressure in hPa (hectopascals).
    """
    # Constants for the Magnus formula
    a = 17.27
    b = 237.7  # in degrees Celsius

    # Magnus formula to calculate saturation vapor pressure
    saturation_vapor_pressure = 6.1078 * math.exp((a * temperature_celsius) / (b + temperature_celsius))

    # Actual vapor pressure
    vapor_pressure = (relative_humidity / 100.0) * saturation_vapor_pressure

    return vapor_pressure


def calculate_dew_point(temperature_celsius: float, vapor_pressure: float) -> float:
    """
    Calculate the dew point temperature given the temperature and vapor pressure.

    Parameters:
    temperature_celsius (float): Temperature in degrees Celsius.
    vapor_pressure (float): Actual vapor pressure in hPa (hectopascals).

    Returns:
    float: Dew point temperature in degrees Celsius.
    """
    # Constants for the Magnus formula
    a = 17.27
    b = 237.7  # in degrees Celsius

    # Calculate saturation vapor pressure at the given temperature
    saturation_vapor_pressure = 6.1078 * math.exp((a * temperature_celsius) / (b + temperature_celsius))

    # Calculate the intermediate value for dew point calculation
    alpha = math.log(vapor_pressure / 6.1078)

    # Dew point calculation
    dew_point = (b * alpha) / (a - alpha)

    return dew_point


@dataclass
class Material:
    name: str
    conductivity: float
    vapour_resistivity: float


@dataclass
class Layer:
    name: str             # Name of the layer
    thickness: float      # Thickness in mm
    material: Material    # Material of the layer

    start_temperature: float = None
    end_temperature: float = None

    start_vapor_pressure: float = None
    end_vapor_pressure: float = None

    def thermal_resistance(self) -> float:
        """
        Calculate the thermal resistance of the layer in m²K/W.
        Converts thickness from mm to meters.
        """
        return (self.thickness / 1000) / self.material.conductivity

    def vapour_resistance(self) -> float:
        """
        Calculate the vapour resistance of the layer in m²sPa/W.
        Converts thickness from mm to meters.
        """
        return (self.thickness / 1000) * self.material.vapour_resistivity


@dataclass
class Wall:
    layers: List[Layer]

    def compute(self, inside_temp: float, outside_temp: float, inside_rh: float, outside_rh: float) -> None:
        total_resistance = sum(layer.thermal_resistance() for layer in self.layers)
        # Calculate the start and end temperatures and vapor pressures
        end_vapor_pressure = calculate_vapor_pressure(outside_temp, outside_rh)
        start_vapor_pressure = calculate_vapor_pressure(inside_temp, inside_rh)
        total_vapor_resistance = sum(layer.vapour_resistance() for layer in self.layers)
        for i, layer in enumerate(self.layers):
            if i == 0:
                layer.start_temperature = inside_temp
                layer.start_vapor_pressure = start_vapor_pressure
            else:
                layer.start_temperature = self.layers[i - 1].end_temperature
                layer.start_vapor_pressure = self.layers[i - 1].end_vapor_pressure

            total_temperature_difference = inside_temp - outside_temp
            this_layer_temperature_difference = total_temperature_difference * (layer.thermal_resistance() / total_resistance)

            layer.end_temperature = layer.start_temperature - this_layer_temperature_difference

            total_vapor_pressure_difference = start_vapor_pressure - end_vapor_pressure
            this_layer_vapor_pressure_difference = total_vapor_pressure_difference * (layer.vapour_resistance() / total_vapor_resistance)
            layer.end_vapor_pressure = layer.start_vapor_pressure - this_layer_vapor_pressure_difference

    def u_value(self) -> float:
        """
        Calculate the U-value of the wall in W/m²K.
        """
        total_resistance = sum(layer.thermal_resistance() for layer in self.layers)

        return 1 / total_resistance

    def display_layers(self) -> None:
        """
        Display the details of each layer.
        """
        for layer in self.layers:
            print(f"{layer.name}: {layer.thickness} mm, {layer.material.conductivity} W/mK (resistance: {layer.thermal_resistance():.4f} m²K/W). Start temperature: {layer.start_temperature:.2f}°C, End temperature: {layer.end_temperature:.2f}°C. Start vapor pressure: {layer.start_vapor_pressure:.2f} hPa, End vapor pressure: {layer.end_vapor_pressure:.2f} hPa.")


# Example layers materials
brick = Material(name="Brickwork", conductivity=0.7, vapour_resistivity=60)
lime_plaster = Material(name="Lime Plaster", conductivity=0.6, vapour_resistivity=5)
rockwool = Material(name="Rockwool", conductivity=0.032, vapour_resistivity=0.5)
wood_wool = Material(name="Wood Wool", conductivity=0.038, vapour_resistivity=5)
pir_insulation = Material(name="PIR Insulation", conductivity=0.022, vapour_resistivity=100)
aluminium_foil = Material(name="Aluminium Foil", conductivity=5, vapour_resistivity=100000)

internal_surface = Layer(name="Internal Surface", thickness=0.01, material=Material(name="Internal Surface", conductivity=(0.01 / (1000 * 0.13)), vapour_resistivity=0))
internal_skim = Layer(name="Internal Lime Plaster Skim", thickness=0, material=lime_plaster)
insulation_board = Layer(name="Insulation Wood Board", thickness=50, material=wood_wool)
internal_render = Layer(name="Internal Lime Plaster Render", thickness=5, material=lime_plaster)
brickwork = Layer(name="Brickwork",  thickness=220, material=brick)
vapor_membrane = Layer(name="Vapor Membrane", thickness=0.4, material=aluminium_foil)
external_surface = Layer(name="External Surface", thickness=0.01, material=Material(name="External Surface", conductivity=(0.01 / (1000 * 0.04)), vapour_resistivity=0))

# Creating a wall with these layers
wall = Wall(layers=[internal_surface, internal_skim, insulation_board, vapor_membrane, internal_render, brickwork, external_surface])

# Calculate and display the U-value
u_value = wall.u_value()
print(f"The U-value of the wall is: {u_value:.3f} W/m²K, {(1.0/u_value):.3f} m²K/W\n\n")  # Industry standard is ~2.1 W/m²K for a solid brick wall

# Optionally display the details of each layer
start_rh = 65
end_rh = 95
outside_temp = 5
inside_temp = 18
wall.compute(inside_temp=inside_temp, outside_temp=outside_temp, inside_rh=start_rh, outside_rh=end_rh)

print(f"Outside temperature: {outside_temp:.1f}°C, Outside relative humidity: {end_rh:.0f}%")
print("Outside dew point temperature: ", calculate_dew_point(outside_temp, calculate_vapor_pressure(outside_temp, end_rh)))

wall.display_layers()

# Plot the temperature profile and the dew point temperature profile
distance = [sum(layer.thickness for layer in wall.layers[:i]) for i in range(len(wall.layers) + 1)]
temperatures = [layer.start_temperature for layer in wall.layers] + [wall.layers[-1].end_temperature]
dew_points = [calculate_dew_point(layer.start_temperature, layer.start_vapor_pressure) for layer in wall.layers] + [calculate_dew_point(wall.layers[-1].end_temperature, wall.layers[-1].end_vapor_pressure)]
vapor_pressures = [layer.start_vapor_pressure for layer in wall.layers] + [wall.layers[-1].end_vapor_pressure]
saturation_vapor_pressures = [calculate_vapor_pressure(layer.start_temperature, 100) for layer in wall.layers] + [calculate_vapor_pressure(wall.layers[-1].end_temperature, 100)]
fig, axs = plt.subplots(1, 2, figsize=(14, 6))
axs[0].plot(distance, temperatures, label="Temperature (°C)")
axs[0].plot(distance, dew_points, label="Dew Point (°C)")
axs[0].set_xlabel("Distance (mm)")
axs[0].set_ylabel("Temperature (°C)")
axs[0].legend()
axs[0].set_title("Condensation Risk Analysis")
axs[0].grid()

axs[1].plot(distance, vapor_pressures, label="Vapor Pressure (hPa)")
axs[1].plot(distance, saturation_vapor_pressures, label="Saturation Vapor Pressure (hPa)")
axs[1].set_xlabel("Distance (mm)")
axs[1].set_ylabel("Vapor Pressure (hPa)")
axs[1].legend()
axs[1].set_title("Vapor Pressure Profile")
axs[1].grid()

# Display the plots
plt.tight_layout()
plt.show()
