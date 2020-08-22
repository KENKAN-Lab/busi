from dataclasses import dataclass
from typing import Union, List


@dataclass
class Side:

    mns: Union[float, int, str, None]
    pls: Union[float, int, str, None]


class SolidLayer:

    def __init__(self, conductivity: float, specific_heat: float, thickness: float, division_number: int):
        """
        Args:
            conductivity: thermal conductivity, W / m K
            specific_heat: volumetric specific heat, J / m3 K
            thickness: thickness, m
            division_number: division_number(int)
        """

        if conductivity <= 0.0:
            raise ValueError('Conductivity shall be more than 0.0.')
        if specific_heat <= 0.0:
            raise ValueError('Specific Heat shall be more than 0.0.')
        if thickness <= 0.0:
            raise ValueError('Thickness shall be more than 0.0.')
        if division_number <= 0.0:
            raise ValueError('Division number shall be more than 0')

        self.__conductivity = conductivity
        self.__specific_heat = specific_heat
        self.__thickness = thickness
        self.__division_number = division_number
        self.__delta_x = thickness / division_number

    @property
    def conductivity(self):
        return self.__conductivity

    @property
    def specific_heat(self):
        return self.__specific_heat

    @property
    def thickness(self):
        return self.__thickness

    @property
    def division_number(self):
        return self.__division_number

    @property
    def delta_x(self):
        return self.__delta_x


class AirLayer:

    def __init__(self, resistance):
        """
        Args:
            resistance: thermal resistance, m2K/W
        """

        if resistance <= 0.0:
            raise ValueError('Air resistance shall be more than 0.0.')

        self.__resistance = resistance

    @property
    def resistance(self):
        return self.__resistance


class Cell:

    def __init__(self, resistance: Side, capacity: float, temp: float):
        """
        Args:
            resistance: Side class
                thermal resistance of side A, m2K/W
                thermal resistance of side B, m2K/W
            capacity: thermal capacity (volumetric specific heat), J/m3K
            temp: temperature, degree C
        """

        self.__resistance: Side = resistance
        self.__capacity: float = capacity
        self.__temp: float = temp

    @property
    def resistance(self):
        return self.__resistance

    @property
    def capacity(self):
        return self.__capacity

    @property
    def temp(self):
        return self.__temp

    @temp.setter
    def temp(self, temp):
        self.__temp = temp

    def get_property(self):
        return {
            'res-': self.__resistance.mns,
            'res+': self.__resistance.pls,
            'cap': self.__capacity,
            'temp': self.__temp,
        }

    @classmethod
    def create(cls, layers: List[Union[SolidLayer, AirLayer]], initial_temp: float):

        if type(layers[0]) == AirLayer:
            raise Exception('First layer should not be air layer.')

        if type(layers[len(layers) - 1]) == AirLayer:
            raise Exception('Last layer should not be air layer.')

        for i, l in enumerate(layers):
            if type(l) == AirLayer:
                if type(layers[i + 1]) == AirLayer:
                    raise Exception('Air layer should not set next to the air layer.')

        cells = []

        for i, layer in enumerate(layers):

            if type(layer) == SolidLayer:

                if i == 0:
                    cells.append(Cell(
                        resistance=Side(
                            mns='surfaceMns',
                            pls=layer.delta_x / layer.conductivity / 2
                        ),
                        capacity=layer.specific_heat * layer.delta_x / 2,
                        temp=initial_temp
                    ))
                elif type(layers[i - 1]) == AirLayer:
                    cells.append(Cell(
                        resistance=Side(
                            mns=layers[i - 1].resistance / 2,
                            pls=layer.delta_x / layer.conductivity / 2
                        ),
                        capacity=layer.specific_heat * layer.delta_x / 2,
                        temp=initial_temp
                    ))
                else:
                    pass

                for j in range(layer.division_number - 1):
                    cells.append(Cell(
                        resistance=Side(
                            mns=layer.delta_x / layer.conductivity / 2,
                            pls=layer.delta_x / layer.conductivity / 2
                        ),
                        capacity=layer.specific_heat * layer.delta_x,
                        temp=initial_temp
                    ))

                if i == len(layers) - 1:
                    cells.append(Cell(
                        resistance=Side(
                            mns=layer.delta_x / layer.conductivity,
                            pls='surfacePls'
                        ),
                        capacity=layer.specific_heat * layer.delta_x / 2,
                        temp=initial_temp
                    ))
                elif type(layers[i + 1]) == AirLayer:
                    cells.append(Cell(
                        resistance=Side(
                            mns=layer.delta_x / layer.conductivity / 2,
                            pls=layers[i + 1].resistance / 2
                        ),
                        capacity=layer.specific_heat * layer.delta_x / 2,
                        temp=initial_temp
                    ))
                else:
                    cells.append(Cell(
                        resistance=Side(
                            mns=layer.delta_x / layer.conductivity / 2,
                            pls=layers[i + 1].delta_x / layers[i + 1].conductivity / 2
                        ),
                        capacity=layer.specific_heat * layer.delta_x / 2
                        + layers[i + 1].specific_heat * layers[i + 1].delta_x / 2,
                        temp=initial_temp
                    ))

        return cells


class Envelope:

    def __init__(self, name: str, cells: List[Cell], long_wave_emissivity: Side):
        """
        Args:
            name: name of Envelope
            cells: cells class
            long_wave_emissivity: Side class
                long wave emissivity of side A
                long wave emissivity of side B
        """

        def check_long_wave_emit_range(value):
            if value < 0.0 or value > 1.0:
                raise ValueError('Long wave emissivity should be in the range between 0.0 and 1.0.')

        check_long_wave_emit_range(long_wave_emissivity.mns)
        check_long_wave_emit_range(long_wave_emissivity.pls)

        self.__name = name
        self.__cells = cells
        self.__e = long_wave_emissivity

    @property
    def cells(self):
        return self.__cells

    @property
    def name(self):
        return self.__name

    @property
    def long_wave_emissivity(self):
        return self.__e

    @property
    def n_cells(self):
        return len(self.__cells)

    @classmethod
    def create(
            cls,
            name: str,
            layers: List[Union[SolidLayer, AirLayer]],
            initial_temp: float = 0.0,
            long_wave_emissivity: Side = Side(mns=0.0, pls=0.0)
    ):
        """
        make Envelope class
        Args:
            name: name
            layers: layers
            initial_temp: initial temperature, degree C
            long_wave_emissivity: Side class
                long wave emissivity of side A
                long wave emissivity of side B
        """

        cells = Cell.create(layers, initial_temp)
        return Envelope(name=name, cells=cells, long_wave_emissivity=long_wave_emissivity)

    def update_temp(self, surf_resistance, air_temp, dt, lwave):

        # temperature of cells, degree C
        temps = [c.temp for c in self.__cells]

        # heat flow minus
        hms = [self._get_heat_flow_mns(air_temp, i, lwave, surf_resistance, temps) for i in range(len(self.__cells))]

        # heat flow plus
        hps = [self._get_heat_flow_pls(air_temp, i, lwave, surf_resistance, temps) for i in range(len(self.__cells))]

        new_temps = [(hm + hp) * dt / c.capacity + t for hm, hp, c, t in zip(hms, hps, self.__cells, temps)]

        for c, t in zip(self.__cells, new_temps):
            c.temp = t

    def _get_heat_flow_mns(self, air_temp, i, lwave, surf_resistance, temps):

        if i == 0:
            return (air_temp.mns - temps[i]) / surf_resistance.mns \
                            + lwave.mns * self.__e.mns - self.get_long_wave(e=self.__e.mns, theta=temps[i])
        else:
            return (temps[i - 1] - temps[i]) / (self.__cells[i - 1].resistance.pls + self.__cells[i].resistance.mns)

    def _get_heat_flow_pls(self, air_temp, i, lwave, surf_resistance, temps):

        if i == len(self.__cells) - 1:
            return (air_temp.pls - temps[i]) / surf_resistance.pls \
                            + lwave.pls * self.__e.pls - self.get_long_wave(e=self.__e.pls, theta=temps[i])
        else:
            return (temps[i + 1] - temps[i]) / (self.__cells[i + 1].resistance.mns + self.__cells[i].resistance.pls)

    @staticmethod
    def get_long_wave(e: float, theta: float) -> float:

        sb_func = 5.670367 * 10 ** (-8)

        tp_water = 273.15

        return e * sb_func * (theta + tp_water) ** 4

    def show_all_temperature(self):
        return [c.temp for c in self.__cells]

    def show_temperature(self, index):
        if index >= len(self.__cells):
            raise IndexError('Index is over the number of cells.')
        if index < 0:
            raise IndexError('Index should be over zero.')

        return self.__cells[index].temp

    def show_surface_temperature(self):

        return Side(self.__cells[0].temp, self.__cells[-1].temp)

    def set_surface_temperature(self, temp: Side):

        if temp.mns is None:
            pass
        else:
            self.__cells[0].temp = temp.mns

        if temp.pls is None:
            pass
        else:
            self.__cells[-1].temp = temp.pls
