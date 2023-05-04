from busi.air.temperature import Temperature
from busi.air import vapor_pressure as vp
from busi.air import saturated_vapor_pressure as svp


class Air:
    """Air class

    Attributes:
        _t(Temperature): temperature, temperature class
        _ah(float): absolute humidity, %
    """

    # The SONNTAG equation is used for calculating the saturated vapour pressure as default.
    svpeq = 'SONNTAG'
    # In this class, as the vapour in the air is treated,
    # 'water' is used as equation type for calculating the saturated vapour pressure as default.
    eqtype = 'water'

    def __init__(self, t, h, humtype):

        self._t = t
        self._ah = None

        self.set_humidity(h, humtype)

    def set_humidity(self, h, humtype):
        """set humidity

        Args:
            h: relative humidity, % or absolute humidity, kg/kgDA
            humtype: 'relative' or 'absolute' as the args type
        """

        f = {
            'relative': self.set_relative_humidity,
            'absolute': self.set_absolute_humidity
        }[humtype]

        f(h)

    # temperature

    def set_temperature(self, t):
        """set temperature

        Args:
            t: temperature, Temperature
        """

        if vp.ah_to_rh(self._ah, t.K, self.svpeq, self.eqtype) > 100.0:
            raise ValueError(
                'Error: Temperature shall be not set at the state that the relative humidity is over 100 %.')

        self._t = t

    def get_temperature(self):
        """get temperature

        Returns: temperature, Temperature
        """

        return self._t

    T = property(get_temperature, set_temperature)

    # relative humidity

    def set_relative_humidity(self, rh: float):
        """set relative humidity

        Args:
            rh: relative humidity, %

        Notes:
            Before this equation is called, the temperature should be set properly.
        """

        if rh < 0.0:
            raise ValueError('Error: Relative humidity shall be equal to or more than 0 %.')

        if rh > 100.0:
            raise ValueError('Error: Relative humidity shall not be over 100 %.')

        self._ah = vp.rh_to_ah(rh, self._t.K, self.svpeq, self.eqtype)

    def get_relative_humidity(self):
        """get relative humidity

        Returns:
            relative humidity, %
        """

        return vp.ah_to_rh(self._ah, self._t.K, self.svpeq, self.eqtype)

    RH = property(get_relative_humidity, set_relative_humidity)

    # absolute humidity

    def set_absolute_humidity(self, ah: float):
        """set absolute humidity

        Args:
            ah: absolute humidity, kg/kgDA
        """

        if ah < 0:
            raise ValueError('Error: Absolute humidity shall be equal to or more than 0 kg/kgDA.')

        if vp.ah_to_rh(ah, self._t.K, self.svpeq, self.eqtype) > 100:
            raise ValueError(
                'Error: Absolute humidity shall be not set at the state that the relative humidity is over 100 %.')

        self._ah = ah

    def get_absolute_humidity(self):
        """get absolute humidity

        Returns:
            absolute humidity, kg/kgDA
        """

        return self._ah

    AH = property(get_absolute_humidity, set_absolute_humidity)

    # vapor pressure

    def set_vapor_pressure(self, pv: float):
        """set vapor pressure

        Args:
            pv: vapor pressure, Pa
        """

        if pv < 0:
            raise ValueError('Error: Vapor pressure shall be equal to or more than 0 Pa.')

        if vp.pv_to_rh(pv, self._t.K, self.svpeq, self.eqtype) > 100:
            raise ValueError(
                'Error: Vapor pressure shall be not set at the state that the relative humidity is over 100 %.')

        self._ah = vp.pv_to_ah(pv)

    def get_vapor_pressure(self):
        """get vapor pressure

        Returns:
            vapour pressure, Pa
        """

        return vp.ah_to_pv(self._ah)

    VP = property(get_vapor_pressure, set_vapor_pressure)

    # saturated vapor pressure

    def get_SVP(self):
        """get the saturated vapor pressure

        Returns:
            saturated vapor pressure, Pa
        """

        return svp.svp(self.svpeq, self.eqtype, self._t.K)

    SVP = property(get_SVP, None)

    # saturated absolute humidity

    def get_SAH(self):
        """get the saturated absolute humidity

        Returns:
            saturated absolute humidiy, kg/kgDA
        """

        return vp._get_saturated_absolute_humidity(self._t.K, self.svpeq, self.eqtype)

    SAH = property(get_SAH, None)

    # dew point temperature

    def get_DPT(self):
        """get the dew point

        Returns:
            dew point, Temperature
        """

        return Temperature(vp.t_dewp(self._ah, self.svpeq, self.eqtype))

    DPT = property(get_DPT, None)
