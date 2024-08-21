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


def calculate_dew_point(vapor_pressure: float) -> float:
    """
    Calculate the dew point temperature given the temperature and vapor pressure.

    Parameters:
    vapor_pressure (float): Actual vapor pressure in hPa (hectopascals).

    Returns:
    float: Dew point temperature in degrees Celsius.
    """
    # Constants for the Magnus formula
    a = 17.27
    b = 237.7  # in degrees Celsius

    # Calculate the intermediate value for dew point calculation
    alpha = math.log(vapor_pressure / 6.1078)

    # Dew point calculation
    dew_point = (b * alpha) / (a - alpha)

    return dew_point


@dataclass
class Material:
    name: str
    conductivity: float
    vapor_resistivity: float


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
        return (self.thickness / 1000) * self.material.vapor_resistivity


# A "parallel" layer: Inherits from Layer, but supports multiple materials each with a different proportion of the total cross section
# usefuly e.g., if you have rafters with rockwool between them
@dataclass
class ParallelLayer:
    name: str
    layers: list[Layer]
    proportions: list[float] # The proportion of each layer, must sum to 1 and have len() == len(layers)

    start_temperature: float = None
    end_temperature: float = None

    start_vapor_pressure: float = None
    end_vapor_pressure: float = None

    @property
    def thickness(self) -> None:
        return self.layers[0].thickness

    def check(self) -> None:
        # Check all thicknesses the same
        assert len(self.layers) > 0
        for layer in self.layers:
            assert layer.thickness == self.layers[0].thickness

        assert sum(self.proportions) == 1.0
        assert len(self.proportions) == len(self.layers)
        # Check all proportions are >= 0.0
        assert all(p >= 0.0 for p in self.proportions)

    def thermal_resistance(self) -> float:
        """
        Calculate the thermal resistance of each layer then combine them
        using the weighted harmonic mean.
        """
        self.check()  # Ensure the proportions and thicknesses are correct
        total_resistance = sum(p / layer.thermal_resistance() for p, layer in zip(self.proportions, self.layers))
        return 1 / total_resistance

    def vapour_resistance(self) -> float:
        """
        Calculate the vapour resistance of each layer then combine them
        using the weighted arithmetic mean.
        """
        self.check()  # Ensure the proportions and thicknesses are correct
        total_resistance = sum(p * layer.vapour_resistance() for p, layer in zip(self.proportions, self.layers))
        return total_resistance


@dataclass
class ExternalConditions:
    inside_temperature: float
    outside_temperature: float

    inside_relative_humidity: float
    outside_relative_humidity: float

    name: str

@dataclass
class Wall:
    name: str
    layers: List[Layer]
    area: float = None  # Area of the wall in m²

    def thickness(self) -> float:
        return sum(layer.thickness for layer in self.layers)

    def compute(self, external_conditions: ExternalConditions) -> None:
        total_resistance = sum(layer.thermal_resistance() for layer in self.layers)
        # Calculate the start and end temperatures and vapor pressures
        end_vapor_pressure = calculate_vapor_pressure(external_conditions.outside_temperature, external_conditions.outside_relative_humidity)
        start_vapor_pressure = calculate_vapor_pressure(external_conditions.inside_temperature, external_conditions.inside_relative_humidity)
        total_vapor_resistance = sum(layer.vapour_resistance() for layer in self.layers)
        for i, layer in enumerate(self.layers):
            if i == 0:
                layer.start_temperature = external_conditions.inside_temperature
                layer.start_vapor_pressure = start_vapor_pressure
            else:
                layer.start_temperature = self.layers[i - 1].end_temperature
                layer.start_vapor_pressure = self.layers[i - 1].end_vapor_pressure

            total_temperature_difference = external_conditions.inside_temperature - external_conditions.outside_temperature
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
            print(f"{layer.name}: {layer.thickness:.2f} mm, {layer.material.conductivity:.2f} W/mK (resistance: {layer.thermal_resistance():.4f} m²K/W). Start temperature: {layer.start_temperature:.2f}°C, End temperature: {layer.end_temperature:.2f}°C. Start vapor pressure: {layer.start_vapor_pressure:.2f} hPa, End vapor pressure: {layer.end_vapor_pressure:.2f} hPa.")

    def display(self, external_conditions: ExternalConditions) -> None:

        # Calculate and display the U-value
        u_value = self.u_value()
        self.compute(external_conditions)

        # Calculate the heat-flux in W
        total_temperature_difference = external_conditions.inside_temperature - external_conditions.outside_temperature
        heat_flux = total_temperature_difference * self.u_value() * self.area

        # Plot the temperature profile and the dew point temperature profile
        distance = [sum(layer.thickness for layer in self.layers[:i]) for i in range(len(self.layers) + 1)]
        temperatures = [layer.start_temperature for layer in self.layers] + [self.layers[-1].end_temperature]
        dew_points = [calculate_dew_point(layer.start_vapor_pressure) for layer in self.layers] + [calculate_dew_point(self.layers[-1].end_vapor_pressure)]
        vapor_pressures = [layer.start_vapor_pressure for layer in self.layers] + [self.layers[-1].end_vapor_pressure]
        saturation_vapor_pressures = [calculate_vapor_pressure(layer.start_temperature, 100) for layer in self.layers] + [calculate_vapor_pressure(self.layers[-1].end_temperature, 100)]

        fig, axs = plt.subplots(1, 2, figsize=(14, 6))
        # Add vertical lines on the plot to indicate the x-position of each layer
        for i in range(len(self.layers)):
            axs[0].axvline(sum(layer.thickness for layer in self.layers[:i]), color="red", linestyle="--", linewidth=0.75)
            axs[1].axvline(sum(layer.thickness for layer in self.layers[:i]), color="red", linestyle="--", linewidth=0.75)

        axs[0].plot(distance, temperatures, label="Temperature (°C)")
        axs[0].plot(distance, dew_points, label="Dew Point (°C)")
        axs[0].set_xlabel("Distance (mm)")
        axs[0].set_ylabel("Temperature (°C)")
        axs[0].legend()
        axs[0].set_title("Condensation Risk Analysis")
        axs[0].grid()

        # If the layer width is more than 20% of the total width, print the layer name at the top of the graph in the mid-x
        # position of the layer
        for i in range(len(self.layers)):
            if self.layers[i].thickness > 0.2 * self.thickness():
                axs[0].text(
                    sum(layer.thickness for layer in self.layers[:i]) + self.layers[i].thickness / 2,
                    max(temperatures + dew_points),
                    self.layers[i].name,
                    ha="center",
                    va="bottom"
                )
                axs[1].text(
                    sum(layer.thickness for layer in self.layers[:i]) + self.layers[i].thickness / 2,
                    max(vapor_pressures + saturation_vapor_pressures),
                    self.layers[i].name,
                    ha="center",
                    va="bottom",
                )

        axs[1].plot(distance, vapor_pressures, label="Vapor Pressure (hPa)")
        axs[1].plot(distance, saturation_vapor_pressures, label="Saturation Vapor Pressure (hPa)")
        axs[1].set_xlabel("Distance (mm)")
        axs[1].set_ylabel("Vapor Pressure (hPa)")
        axs[1].legend()
        axs[1].set_title("Vapor Pressure Profile")
        axs[1].grid()
        fig.suptitle(f"Temperature and Vapor Pressure Profile for '{self.name}' with a U-value of {u_value:.3f} W/m²K in conditions: {external_conditions.name}. Heat Flux: {heat_flux:,.0f} W.")

        # Add text in the top right to say "Outside"
        axs[1].text(
            sum(layer.thickness for layer in self.layers[:i]) + self.layers[i].thickness / 2,
            max(vapor_pressures + saturation_vapor_pressures) + 0.05,
            "Outside",
            ha="center",
            va="top",
        )

        # Add text in the top left to say "Inside"
        axs[1].text(
            0,
            max(vapor_pressures + saturation_vapor_pressures) + 0.05,
            "Inside",
            ha="center",
            va="top",
        )

        # Display the plots
        plt.tight_layout()
        plt.show()


# Example layers materials
brick = Material(name="Brickwork", conductivity=0.8, vapor_resistivity=60)
lime_plaster = Material(name="Lime Plaster", conductivity=0.5, vapor_resistivity=5)
rockwool = Material(name="Rockwool", conductivity=0.038, vapor_resistivity=0.5)
wood_wool = Material(name="Wood Wool", conductivity=0.038, vapor_resistivity=5)
pir_insulation = Material(name="PIR Insulation", conductivity=0.022, vapor_resistivity=100)
phenolic_insulation = Material(name="Phenolic Insulation", conductivity=0.019, vapor_resistivity=400)
metalized_foil = Material(name="Aluminium Foil", conductivity=5, vapor_resistivity=100000)
plasterboard_material = Material(name="Plasterboard", conductivity=0.25, vapor_resistivity=50)
wood = Material(name="Wood", conductivity=0.13, vapor_resistivity=50)
osb = Material(name="OSB", conductivity=0.13, vapor_resistivity=50)

internal_surface = Layer(name="Internal Surface", thickness=0.01, material=Material(name="Internal Surface", conductivity=(0.01 / (1000 * 0.13)), vapor_resistivity=0))
internal_skim = Layer(name="Internal Lime Plaster Skim", thickness=5, material=lime_plaster)
insulation_board = Layer(name="Insulation Wood Board", thickness=100, material=wood_wool)
insulation_board_40 = Layer(name="Insulation Wood Board, 40mm", thickness=40, material=wood_wool)

internal_render = Layer(name="Internal Lime Plaster Render", thickness=3, material=lime_plaster)
internal_original_plaster = Layer(name="Internal Original Plaster", thickness=18, material=lime_plaster)
brickwork = Layer(name="Brickwork",  thickness=220, material=brick)
vapor_membrane = Layer(name="Vapor Membrane", thickness=0.2, material=metalized_foil)
external_surface = Layer(name="External Surface", thickness=0.01, material=Material(name="External Surface", conductivity=(0.01 / (1000 * 0.04)), vapor_resistivity=0))

pir_insulation_40 = Layer(name="PIR Insulation, 40mm", thickness=40, material=pir_insulation)
pir_insulation_50 = Layer(name="PIR Insulation, 50mm", thickness=50, material=pir_insulation)
pir_insulation_25 = Layer(name="PIR Insulation, 25mm", thickness=25, material=pir_insulation)
phenolic_insulation_40 = Layer(name="Phenolic Insulation, 40mm", thickness=40, material=phenolic_insulation)
plasterboard = Layer(name="Plasterboard", thickness=12.5, material=plasterboard_material)

# Creating a wall with these layers
wall_original = Wall(name="Original Wall", layers=[internal_surface, internal_original_plaster, brickwork, external_surface], area=6*2.5)
wall_40 = Wall(name="Steico 40/Lime", layers=[internal_surface, internal_skim, insulation_board_40, internal_render, brickwork, external_surface], area=6*2.5)
wall_100 = Wall(name="Steico 100/Lime", layers=[internal_surface, internal_skim, insulation_board, internal_render, brickwork, external_surface], area=6*2.5)
wall_2 = Wall(name="Steico 100/Lime w/Vapor Barrier", layers=[internal_surface, internal_skim, vapor_membrane, insulation_board, internal_render, brickwork, external_surface], area=6*2.5)
wall_pir = Wall(name="PIR Insulation", layers=[internal_surface, plasterboard, pir_insulation_40, brickwork, external_surface], area=6*2.5)
wall_pir_25 = Wall(name="PIR Insulation 25mm", layers=[internal_surface, plasterboard, pir_insulation_25, internal_render, brickwork, external_surface], area=6*2.5)
wall_pir_50 = Wall(name="PIR Insulation 50mm", layers=[internal_surface, plasterboard, pir_insulation_50, internal_render, brickwork, external_surface], area=6*2.5)
wall_phenolic = Wall(name="Phenolic Insulation", layers=[internal_surface, plasterboard, phenolic_insulation_40, internal_render, brickwork, external_surface], area=6*2.5)

# Typical internal/external conditions for UK
summer_day = ExternalConditions(inside_temperature=18, outside_temperature=25, inside_relative_humidity=65, outside_relative_humidity=65, name="Summer Day")
summer_night = ExternalConditions(inside_temperature=17, outside_temperature=12, inside_relative_humidity=65, outside_relative_humidity=65, name="Summer Night")
wet_winter_day = ExternalConditions(inside_temperature=18, outside_temperature=4, inside_relative_humidity=65, outside_relative_humidity=95, name="Wet Winter Day")
dry_winter_day = ExternalConditions(inside_temperature=18, outside_temperature=-3, inside_relative_humidity=65, outside_relative_humidity=55, name="Dry Winter Day")

wall_pir_50.compute(wet_winter_day)
wall_pir_50.display_layers()

wall_original.display(wet_winter_day)
wall_40.display(wet_winter_day)
# wall_100.display(wet_winter_day)
wall_pir.display(wet_winter_day)
wall_pir_25.display(wet_winter_day)
wall_phenolic.display(wet_winter_day)
wall_pir_50.display(wet_winter_day)

# Using the U value (W/m2K) calculate the temperature gradient for a given amount of heat input (200W)
power_input = 100
wall_original_dt = power_input / (wall_original.u_value() * wall_original.area)
wall_40_dt = power_input / (wall_40.u_value() * wall_40.area)
wall_pir_dt = power_input / (wall_pir.u_value() * wall_pir.area)

print(f"Temperature gradient for Original Wall: {wall_original_dt:.2f}°C")
print(f"Temperature gradient for Steico 40/Lime: {wall_40_dt:.2f}°C")
print(f"Temperature gradient for PIR Insulation: {wall_pir_dt:.2f}°C")


original_ceiling = Wall(
    name="Original ceiling",
    layers=[internal_surface, internal_original_plaster, external_surface],
    area=6*4
)

rockwool_200 = ParallelLayer(
    name="Rockwool to Joists",
    layers=[
        Layer(name="Rockwool", thickness=200, material=rockwool),
        Layer(name="Joists", thickness=200, material=wood)
    ],
    proportions=[
        0.9,
        0.1
    ]
)
rockwool_200_ceiling = Wall(
    name="Rockwool to Joists",
    layers=[
        internal_surface,
        internal_original_plaster,
        rockwool_200,
        external_surface
    ],
    area=6*4
)
rockwool_300_ceiling = Wall(
    name="Gov recommended rockwool",
    layers=[
        internal_surface,
        internal_original_plaster,
        rockwool_200,
        Layer(name="Rockwool", thickness=100, material=rockwool),
        external_surface
    ],
    area=6*4
)

rockwool_plasterboard = Wall(
    name="Rockwool & insulated plasterboard",
    layers=[
        internal_surface,
        plasterboard,
        pir_insulation_50,
        rockwool_200,
        external_surface
    ],
    area=6*4
)

original_ceiling.display(wet_winter_day)
rockwool_200_ceiling.display(wet_winter_day)
rockwool_300_ceiling.display(wet_winter_day)
rockwool_plasterboard.display(wet_winter_day)
