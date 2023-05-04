class Temperature:
    """
    keep the temperature.
    Temperature is represented by the two units, one is 'Kelvin' and the other is 'Celsius'.
    This class can be convert the units automatically when getting or setting by specifying the unit type.

    Attributes:
        kelvin(float): temperature, which unit is Kelvin, K
    """

    def __init__(self, v, unit='Kelvin'):
        """
        Args:
            v: value
            unit: 'Kelvin' or 'Celsius'
        """
        self.kelvin = None
        self.set_value(v, unit)

    def set_as_celsius(self, v):
        """set the value as celsius

        Args:
            v: temperature, degree C
        """
        self.set_value(v, 'Celsius')

    def get_as_celsius(self):
        """get the value as celsius

        Returns:
            temperature, degree C
        """

        return self.get_value('Celsius')

    def set_as_kelvin(self, v):
        """set the value as kelvin

        Args:
            v: temperature, K
        """
        self.set_value(v, 'Kelvin')

    def get_as_kelvin(self):
        """get the value as kelvin

        Returns:
            temperature, K
        """
        return self.get_value('Kelvin')

    def set_value(self, v, unit):
        """set value as Celsius or Kelvin

        Args:
            v: value
            unit: 'Kelvin' or 'Celsius'
        """

        self.kelvin = {
            'Kelvin': v,
            'Celsius': v + 273.15
        }[unit]

        if self.kelvin < 0:
            raise ValueError('Error: Temperature should be equal to or more than 0 kelvin.')

    def get_value(self, unit):
        """get value as Celsius or Kelvin

        Args:
            unit: the type of the units

        Returns:
            temperature, K or degree C
        """

        return {
            'Kelvin': self.kelvin,
            'Celsius': self.kelvin - 273.15
        }[unit]

    # set property
    C = property(get_as_celsius, set_as_celsius)
    K = property(get_as_kelvin, set_as_kelvin)

    def __add__(v1, v2):
        return Temperature(v1.K + v2)

    def __sub__(v1, v2):
        return Temperature(v1.K - v2)

    def difference(v1, v2):
        return v1.K - v2.K

    def __lt__(v1, v2):
        return v1.K < v2.K

    def __le__(v1, v2):
        return v1.K <= v2.K

    def __eq__(v1, v2):
        return v1.K == v2.K

    def __ne__(v1, v2):
        return v1.K != v2.K

    def __ge__(v1, v2):
        return v1.K >= v2.K

    def __gt__(v1, v2):
        return v1.K > v2.K
