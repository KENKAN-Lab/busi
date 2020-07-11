""" Evaluation on thermal resistance of float glasses and thermal transmittance of glazing

This module is for evaluation on thermal resistance of float glassess and thermal transmittance of glazing, which is
based on JIS R3107 1988.

"""

from enum import Enum


class GlassDirection(Enum):
    """Direction of the glass

    """
    # Air layer is vertical and heat flow direction is horizontal.
    VERTICAL = 1
    # Air layer is horizontal and heat flow direction is upward.
    HORIZONTAL = 2
    # Air layer faces 45 degree and heat flow direction is upward.
    SLOPE = 3


class AirType(Enum):
    """Type of air

    """

    AIR = 1
    ARGON = 2
    SF6 = 3
    KRYPTON = 4

# region air property functions


def get_interporation_property(temp, v1, v2, v3, v4):
    """
    define the physical property at the designated temperature
        Args:
            temp: temperature, degree C
            v1: the physical property at -10 degree C
            v2: the physical property at 0 degree C
            v3: the physical property at 10 degree C
            v4: the physical property at 20 degree C

        Returns:
            the physical property
    """

    t1 = -10.0
    t2 = 0.0
    t3 = 10.0
    t4 = 20.0

    if temp < t1:
        raise ValueError('空気層の指定は-10℃から20℃の間でなければなりません')
    elif temp < t2:
        return v1 * (temp - t2) / (t1 - t2) + v2 * (temp - t1) / (t2 - t1)
    elif temp < t3:
        return v2 * (temp - t3) / (t2 - t3) + v3 * (temp - t2) / (t3 - t2)
    elif temp < t4:
        return v3 * (temp - t4) / (t3 - t4) + v4 * (temp - t3) / (t4 - t3)
    else:
        raise ValueError('空気層の指定は-10℃から20℃の間でなければなりません')


def get_rho(temp: float, air_type: AirType) -> float:
    """
    calculate the density of the air

    Args:
        temp: temperature, degree C
        air_type: type of the air

    Returns:
        density, kg/m3
    """

    rho1, rho2, rho3, rho4 = {
        AirType.AIR: (1.326, 1.277, 1.232, 1.189),
        AirType.ARGON: (1.829, 1.762, 1.699, 1.640),
        AirType.SF6: (6.844, 6.602, 6.360, 6.118),
        AirType.KRYPTON: (3.832, 3.690, 3.560, 3.430)
    }[air_type]

    return get_interporation_property(temp, rho1, rho2, rho3, rho4)


def get_mu(temp: float, air_type: AirType) -> float:
    """
    calculate the viscosity(mu value) of the air
    Args:
        temp: temperature, degree C
        air_type: type of the air

    Returns:
        viscosity(mu value), kg/ms
    """

    mu1, mu2, mu3, mu4 = {
        AirType.AIR: (1.661, 1.711, 1.761, 1.811),
        AirType.ARGON: (2.038, 2.101, 2.164, 2.228),
        AirType.SF6: (1.383, 1.421, 1.459, 1.497),
        AirType.KRYPTON: (2.260, 2.330, 2.400, 2.470)
    }[air_type]

    return get_interporation_property(temp, mu1, mu2, mu3, mu4) * 10 ** -5


def get_lambda(temp: float, air_type: AirType) -> float:
    """
    calculate thermal conductivity of the air
    Args:
        temp: temperature, degree C
        air_type: type of the air

    Returns:
        thermal conductivity, W/mK
    """

    lambda1, lambda2, lambda3, lambda4 = {
        AirType.AIR: (2.336, 2.416, 2.496, 2.576),
        AirType.ARGON: (1.584, 1.634, 1.684, 1.734),
        AirType.SF6: (1.119, 1.197, 1.275, 1.354),
        AirType.KRYPTON: (0.842, 0.870, 0.900, 0.926)
    }[air_type]

    return get_interporation_property(temp, lambda1, lambda2, lambda3, lambda4) * 10 ** -2


def get_c_air(air_type: AirType) -> float:
    """
    calculate the specific heat of the air

    Args:
        air_type: type of the air

    Returns:
        specific heat, J/kgK
    """

    return {
        AirType.AIR: 1.008 * 10 ** 3,
        AirType.ARGON: 0.519 * 10 ** 3,
        AirType.SF6: 0.614 * 10 ** 3,
        AirType.KRYPTON: 0.245 * 10 ** 3
    }[air_type]

# endregion


ABSOLUTE_ZERO = 273.15


def get_sigma():
    """
    get Stefan Boltzmann coefficient

    Returns:
        Stefan Boltzmann coefficient
    """

    return 5.67 / 100000000


def get_cc(nu: float) -> float:
    """
    calculate convective effect coefficient

    Args:
        nu: Nusselt number

    Returns:
        convective effect coefficient

    Notes:
        equation (5)
    """

    return min(1.0, nu)


def get_nu(a: float, n: float, gr: float, pr: float) -> float:
    """
    calculate Nusselt value

    Args:
        a: a value
        n: n value
        gr: Grasshof number
        pr: Prandtl number

    Returns:
        Nusselt value

    Notes:
        equation(6)
    """

    return a * (gr * pr)**n


def get_gr(delta_t: float, t_dash_m: float, s: float, rho: float, mu_air) -> float:
    """get the Grasshof number

    Args:
        delta_t: temperature difference of the surfaces of the air layer, K
        t_dash_m: average temperature of the air of the air layer, K
        s: thickness, m
        rho: density, kg/m3
        mu_air: viscosity(mu value), kg/ms

    Returns:
        Grasshof number

    Notes:
        equation (7)
    """

    return 9.81 * s ** 3 * delta_t * rho ** 2 / t_dash_m / mu_air ** 2


def get_pr(mu_air: float, c_air: float, lambda_air: float) -> float:
    """
    get Prandtl number

    Args:
        mu_air: viscosity(mu value), kg/ms
        c_air: specific heat, J/kgK
        lambda_air: thermal conductivity, W/mK

    Returns:
        Prandtl number

    Notes:
        equation (8)
    """

    return mu_air * c_air / lambda_air


def get_a(direction: GlassDirection) -> float:
    """

    Returns:
        value A

    Notes:
        5.3.1
    """

    return {
        GlassDirection.VERTICAL: 0.035,
        GlassDirection.HORIZONTAL: 0.16,
        GlassDirection.SLOPE: 0.10
    }[direction]


def get_n(direction: GlassDirection) -> float:
    """

    Returns:
        value n

    Notes:
        5.3.1
    """

    return {
        GlassDirection.VERTICAL: 0.38,
        GlassDirection.HORIZONTAL: 0.28,
        GlassDirection.SLOPE: 0.31
    }[direction]


def get_h_g(direction, delta_t, t_dash_m, s, rho_air, mu_air, c_air, lambda_air):
    """
    get air conductance
    Args:
        direction: direction
        delta_t: temperature difference of the surfaces of the air layer, K
        t_dash_m: average temperature of the air of the air layer, K
        s:
        rho_air:
        mu_air:
        c_air:
        lambda_air:

    Returns:
        air conductance, W/m2K

    Notes:
        equation (5)
    """

    # a value
    a = get_a(direction)

    # n value
    n = get_n(direction)

    # Grasshof number
    gr = get_gr(delta_t, t_dash_m, s, rho_air, mu_air)

    # Prandtl number
    pr = get_pr(mu_air, c_air, lambda_air)

    # Nusselt value
    nu = get_nu(a, n, gr, pr)

    # convective effect coefficient
    cc = get_cc(nu)

    return cc * lambda_air / s


def get_h_r(t_m: float, emissivity1: float, emissivity2: float) -> float:
    """
    calculate radiative conductance
    Args:
        t_m: average of the absolute temperatures at each glass surface of the air layer, K
        emissivity1: emissivity of glass surface 1
        emissivity2: emissivity of glass surface 2

    Returns:
        radiative conductance, W/m2K

    Notes:
        equation (4)
    """

    return 4 * get_sigma() / (1 / emissivity1 + 1 / emissivity2 - 1) * t_m**3


class Surface:
    """Surface on each side of the air layer

    """

    def __init__(self):
        pass

    @property
    def emissivity(self):
        """get the corrected emissivity

        Returns:
            corrected emissivity

        Notes:
            This method is abstract method and it should be overrided..
        """
        raise NotImplementedError()


class NonLowESurface(Surface):
    """Surface of the glass with non low-E thin film

    """

    def __init__(self):
        super().__init__()

    @property
    def emissivity(self):
        """get the corrected emissivity

        Returns:
            corrected emissivity

        Notes:
            section 5.2
        """

        return 0.837


class LowESurface(Surface):
    """Surface of the glass with low-E thin film

    """

    def __init__(self, epsilon_n: float):
        """

        Args:
            epsilon_n: vertical emissivity
        """

        self._epsilon_n = epsilon_n
        super().__init__()

    @property
    def epsilon_n(self):

        return self._epsilon_n

    @property
    def emissivity(self):
        """get the corrected emissivity

        Returns:
            corrected emissivity

        Notes:
            section 5.2
        """

        x1, y1 = 0.03, 1.22
        x2, y2 = 0.05, 1.18
        x3, y3 = 0.1, 1.14
        x4, y4 = 0.2, 1.1
        x5, y5 = 0.3, 1.06
        x6, y6 = 0.4, 1.03
        x7, y7 = 0.5, 1.0
        x8, y8 = 0.6, 0.98
        x9, y9 = 0.7, 0.96
        x10, y10 = 0.8, 0.95
        x11, y11 = 0.89, 0.94

        # vertical emissivity
        e = self.epsilon_n

        if e < x1:
            raise ValueError()
        elif e < x2:
            return y1 * (x2 - e) / (x2 - x1) + y2 * (e - x1) / (x2 - x1)
        elif e < x3:
            return y2 * (x3 - e) / (x3 - x2) + y3 * (e - x2) / (x3 - x2)
        elif e < x4:
            return y3 * (x4 - e) / (x4 - x3) + y4 * (e - x3) / (x4 - x3)
        elif e < x5:
            return y4 * (x5 - e) / (x5 - x4) + y5 * (e - x4) / (x5 - x4)
        elif e < x6:
            return y5 * (x6 - e) / (x6 - x5) + y6 * (e - x5) / (x6 - x5)
        elif e < x7:
            return y6 * (x7 - e) / (x7 - x6) + y7 * (e - x6) / (x7 - x6)
        elif e < x8:
            return y7 * (x8 - e) / (x8 - x7) + y8 * (e - x7) / (x8 - x7)
        elif e < x9:
            return y8 * (x9 - e) / (x9 - x8) + y9 * (e - x8) / (x9 - x8)
        elif e < x10:
            return y9 * (x10 - e) / (x10 - x9) + y10 * (e - x9) / (x10 - x9)
        elif e < x11:
            return y10 * (x11 - e) / (x11 - x10) + y11 * (e - x10) / (x11 - x10)
        else:
            raise ValueError()


class AirProperty:
    """Air in the air layer

    """

    def __init__(self, c_air: float, c_argon: float, c_sf6: float, c_krypton):
        """

        Args:
            c_air: volumetric ratio of air, %
            c_argon: volumetric ratio of argon, %
            c_sf6: volumetric ratio of SF6, %
            c_krypton: volumetric ratio of krypton, %
        """

        # confirm that the sum of the volumetric ratio of each air property equals to 100%.
        if c_air + c_argon + c_sf6 + c_krypton != 100.0:
            raise ValueError('気体の容積割合の合計は100.0になる必要があります')

        self._c_air = c_air
        self._c_argon = c_argon
        self._c_sf6 = c_sf6
        self._c_krypton = c_krypton

    def get_mixed_property(self, f):

        return f(AirType.AIR) * self._c_air\
            + f(AirType.ARGON) * self._c_argon\
            + f(AirType.SF6) * self._c_sf6\
            + f(AirType.KRYPTON) * self._c_krypton

    def get_rho(self, temp: float) -> float:
        """
        calculate the density of the air
        Args:
            temp: temperature, degree C

        Returns:
            density, kg/m3
        """

        return self.get_mixed_property(lambda air_type: get_rho(temp, air_type))

    def get_mu(self, temp: float) -> float:
        """
        calculate the viscosity(mu value) of the air
        Args:
            temp: temperature, degree C

        Returns:
            viscosity(mu value), kg/ms
        """

        return self.get_mixed_property(lambda air_type: get_mu(temp, air_type))

    def get_lambda(self, temp: float) -> float:
        """
        calculate thermal conductivity of the air
        Args:
            temp: temperature, degree C

        Returns:
            thermal conductivity, W/mK
        """

        return self.get_mixed_property(lambda air_type: get_lambda(temp, air_type))

    def get_c(self, temp: float) -> float:
        """
        calculate the specific heat of the air
        Args:
            temp: temperature, degree C

        Returns:
            specific heat, J/kgK
        """

        return self.get_mixed_property(lambda air_type: get_c_air(air_type))


class AirLayer:

    def __init__(
            self,
            surface1: Surface, surface2: Surface,
            air_property: AirProperty,
            direction: GlassDirection,
            s: float
    ):
        """

        Args:
            surface1: glass surface 1
            surface2: glass surface 2
            air_property: air property
            direction: direction of the air layer
            s: thickness, m
        """

        self._surface1 = surface1
        self._surface2 = surface2
        self._air_property = air_property
        self._direction = direction
        self._s = s

    def get_hs(self, t_m: float, delta_t: float, t_dash_m: float) -> float:
        """
        calculate thermal conductance
        Args:
            t_m: average of the absolute temperatures at each glass surface of the air layer, K
            delta_t: temperature difference of the surfaces of the air layer, K
            t_dash_m: average temperature of the air of the air layer, K

        Returns:
            thermal conductance, W/m2K

        Notes:
            equation (2)
        """

        return self.get_h_r(t_m) + self.get_h_g(delta_t, t_dash_m)

    def get_h_r(self, t_m: float) -> float:
        """
        calculate radiative conductance
        Args:
            t_m: average of the absolute temperatures at each glass surface of the air layer, K

        Returns:
            radiative conductance, W/m2K

        Notes:
            equation (4)
        """

        return get_h_r(t_m=t_m, emissivity1=self._surface1.emissivity, emissivity2=self._surface2.emissivity)

    def get_h_g(self, delta_t: float, t_dash_m: float) -> float:
        """
        calculate the air conductance
        Args:
            delta_t: temperature difference of the surfaces of the air layer, K
            t_dash_m: average temperature of the air of the air layer, K

        Returns:
            air conductance, W/m2K
        """

        theta_dash_m = t_dash_m - ABSOLUTE_ZERO

        return get_h_g(
            direction=self._direction,
            delta_t=delta_t,
            t_dash_m=t_dash_m,
            s=self._s,
            rho_air=self._air_property.get_rho(theta_dash_m),
            mu_air=self._air_property.get_mu(theta_dash_m),
            c_air=self._air_property.get_c(theta_dash_m),
            lambda_air=self._air_property.get_lambda(theta_dash_m))


class GlassLayer:
    """
    This class represents the unit glass.

    Attributes:
        _t_g: thickness of the glass, m
        _lambda_g: thermal conductivity of the glass, W/mK (= 1.0)
    """

    def __init__(self, t_g: float):
        """
        Args:
            t_g: thickness of glass, m
        """

        self._t_g = t_g
        self._lambda_g = 1.0

    @property
    def t_g(self):
        """
        Returns:
            thickness of glass, m
        """

        return self._t_g

    @property
    def lambda_g(self):
        """
        Returns:
            thermal conductivity, W/mK
        Notes:
            section 5.1
        """

        return self._lambda_g


class Glass:

    def __init__(self):
        pass

    @staticmethod
    def get_epsilon():

        # section 5.2
        # This value should be used only for the surface of the inside and outside, which is not the low-e surface.
        return 0.837

    @staticmethod
    def get_h_i():
        """

        Returns:
            heat transfer coefficient on the inside surface, W/m2 K
        Notes:
            This value is supposed to apply for the vertical glasses.
        """

        # section 5.4.2
        return 5.4 * Glass.get_epsilon() + 4.1

    @staticmethod
    def get_h_e():
        """

        Returns:
            heat transfer coefficient on the outside surface, W/m2 K
        Notes:
            This value is supposed to apply for the vertical glasses.
        """

        # section 5.4.2
        return 4.9 * Glass.get_epsilon() + 16.3

    def get_heat_transmittance(self):
        """
        calculate the heat transmittance

        Returns:
            heat transmittance, W/m2K

        Notes:
            equation (3)
        """

        return 1 / (1 / Glass.get_h_i() + self.get_heat_resistance() + 1 / Glass.get_h_e())

    def get_heat_resistance(self) -> float:
        """
            get the heat resistance
            This method is override at the upper class.

        Returns:
            heat resistance, m2K/W

        Notes:
            This equation should be override.
        """

        raise NotImplementedError()


class GlassSingle(Glass):

    def __init__(self, gl: GlassLayer):

        self._gl = gl
        super().__init__()

    def get_heat_resistance(self):
        """get heat resistance

        Returns:
            heat resistance of glass, m2K/W

        Notes:
            equation (1)
        """

        return self._gl.t_g / self._gl.lambda_g


class GlassDouble(Glass):

    def __init__(self, gl1: GlassLayer, al: AirLayer, gl2: GlassLayer):
        self._gl1 = gl1
        self._gl2 = gl2
        self._al = al
        super().__init__()

    def get_heat_resistance(self):
        """
        get heat resistance
        Returns:
            heat resistance, m2K/W

        Notes:
            equation (1)
        """

        # section 5.4.1
        # average of the absolute temperatures at each glass surface of the air layer, K
        t_m = 283.0

        # temperature difference of the surfaces of the air layer, K
        delta_t = 15.0

        # average temperature of the air of the air layer, K
        t_dash_m = 283.0

        return self._gl1.t_g / self._gl1.lambda_g + self._al.get_hs(t_m, delta_t, t_dash_m) + self._gl2.t_g / self._gl2.lambda_g


class GlassTriple(Glass):

    def __init__(self, gl1: GlassLayer, al1: AirLayer, gl2: GlassLayer, al2: AirLayer, gl3: GlassLayer):
        self._gl1 = gl1
        self._al1 = al1
        self._gl2 = gl2
        self._al2 = al2
        self._gl3 = gl3
        super().__init__()

    def get_heat_resistance(self):
        """
        get heat resistance
        Returns:
            heat resistance, m2K/W

        Notes:
            equation (1)
        """

        # section 5.4.1
        # average of the absolute temperatures at each glass surface of the air layer, K
        t_m = 283.0

        # temperature difference of the surfaces of the air layer, K
        delta_t = 15.0

        # average temperature of the air of the air layer, K
        t_dash_m = 283.0

        return self._gl1.t_g / self._gl1.lambda_g\
            + self._al1.get_hs(t_m, delta_t, t_dash_m)\
            + self._gl2.t_g / self._gl2.lambda_g\
            + self._al2.get_hs(t_m, delta_t, t_dash_m)\
            + self._gl3.t_g / self._gl3.lambda_g


if __name__ == '__main__':

    g1 = GlassSingle(gl=GlassLayer(t_g=0.003))
    u1 = g1.get_heat_transmittance()
    print(u1)

    g2 = GlassDouble(
        gl1=GlassLayer(t_g=0.003),
        al=AirLayer(
            surface1=NonLowESurface(),
            surface2=NonLowESurface(),
            air_property=AirProperty(c_air=100.0, c_argon=0.0,c_sf6=0.0,c_krypton=0.0),
            direction=GlassDirection.VERTICAL,
            s=0.007),
        gl2=GlassLayer(t_g=0.003)
    )
    u2 = g2.get_heat_transmittance()
    print(u2)

    g3 = GlassTriple(
        gl1=GlassLayer(t_g=0.003),
        al1=AirLayer(
            surface1=NonLowESurface(),
            surface2=NonLowESurface(),
            air_property=AirProperty(c_air=100.0, c_argon=0.0,c_sf6=0.0,c_krypton=0.0),
            direction=GlassDirection.VERTICAL,
            s=0.007),
        gl2=GlassLayer(t_g=0.003),
        al2 = AirLayer(
            surface1=NonLowESurface(),
            surface2=NonLowESurface(),
            air_property=AirProperty(c_air=100.0, c_argon=0.0, c_sf6=0.0, c_krypton=0.0),
            direction=GlassDirection.VERTICAL,
            s=0.007),
        gl3 = GlassLayer(t_g=0.003)
    )
    u3 = g3.get_heat_transmittance()
    print(u3)