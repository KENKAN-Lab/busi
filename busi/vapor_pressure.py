import scipy.optimize

import src.saturated_vapor_pressure as svp


def get_atmospheric_pressure() -> float:
    """get atmospheric pressure

    Returns:
        atmospheric pressure, Pa
    """

    return 101325.0


def ah_to_pv(x: float) -> float:
    """calculate the vapour pressure from the absolute humidity

    Args:
        x: absolute humidity, kg/kgDA

    Returns:
        vapour pressure, Pa
    """

    p_atm = get_atmospheric_pressure()

    return p_atm * x / (0.622 + x)


def pv_to_ah(p_v: float) -> float:
    """calculate the absolute humidity from the vapour pressure

    Args:
        p_v: vapour pressure, Pa

    Returns:
        absolute humidity, kg/kgDA
    """

    p_atm = get_atmospheric_pressure()

    return 0.622 * p_v / (p_atm - p_v)


def pv_to_rh(p_v: float, t: float, equation: str, status: str) -> float:
    """calculate the relative humidity from the vapour pressure and temperature

    Args:
        p_v: vapour pressure, Pa
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        relative humidity, %
    """

    return p_v / svp.get_saturated_vapor_pressure(equation=equation, status=status, t=t) * 100


def rh_to_pv(phi: float, t: float, equation: str, status: str) -> float:
    """calculate the vapour pressure from the relative humidity and temperature

    Args:
        phi: relative humidity, %
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        vapour pressure, Pa
    """

    return svp.get_saturated_vapor_pressure(equation=equation, status=status, t=t) * phi / 100


def rh_to_ah(phi: float, t: float, equation: str, status: str) -> float:
    """calculate absolute humidity from relative humidity and temperature

    Args:
        phi: relative humidity, %
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        absolute humidity, kg/kgDA
    """

    p_v = rh_to_pv(phi, t, equation, status)

    return pv_to_ah(p_v)


def ah_to_rh(x: float, t: float, equation: str, status: str) -> float:
    """calculate relative humidity from absolute humidity and temperature

    Args:
        x: absolute humidity, kg/kgDA
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        relative humidity, %
    """

    p_v = ah_to_pv(x)

    return pv_to_rh(p_v, t, equation, status)


def get_saturated_absolute_humidity(t: float, equation: str, status: str) -> float:
    """calculate the saturated absolute humidity from the temperature

    Args:
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        saturated absolute humidity, kg/kgDA
    """

    p_sv = svp.get_saturated_vapor_pressure(equation, status, t)

    return pv_to_ah(p_sv)


def get_dew_point_temp(x_s: float, equation: str, status: str) -> float:
    """calculate the dew point temperature from the absolute humidity

    Args:
        x_s: absolute humidity, kg/kgDA
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        dew point temperature, K
    """

    def f(t):
        return get_saturated_absolute_humidity(t, equation, status) - x_s

    return scipy.optimize.newton(f, 273.15)

    # use Bisection method
#    TL = 0.0    # T = 0.0 (K)
#    TH = 373.15 # T = 373.15 (K) = 100.0 (degree C)
#    n = 0 # Counter
#    nmax = 10000 # Maximum number of attempts
#    while True:
#        T = ( TH + TL )/2
#        sah = get_saturated_absolute_humidity(T, equation, status)
#        if abs(sah - x_s) < 1.0e-7:
#            break
#        n = n + 1
#        if n >= nmax:
#            raise Error('Error: Not converged in the calculation of the dew point temperature.')
#        if x_s < sah:
#            TH = T
#        else:
#            TL = T
#    return T
