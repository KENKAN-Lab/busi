import scipy.optimize

import busi.air.saturated_vapor_pressure as svp


def get_atmospheric_pressure() -> float:
    """Get atmospheric pressure.

    Returns:
        atmospheric pressure, Pa
    """

    return 101325.0


def ah_to_pv(x: float) -> float:
    """Calculate the vapour pressure from the absolute humidity.

    Args:
        x: absolute humidity, kg/kgDA

    Returns:
        vapour pressure, Pa
    """

    p_atm = get_atmospheric_pressure()

    return p_atm * x / (0.622 + x)


def pv_to_ah(p_v: float) -> float:
    """Calculate the absolute humidity from the vapour pressure.

    Args:
        p_v: vapour pressure, Pa

    Returns:
        absolute humidity, kg/kgDA
    """

    p_atm = get_atmospheric_pressure()

    return 0.622 * p_v / (p_atm - p_v)


def pv_to_rh(p_v: float, t: float, equation: str, status: str) -> float:
    """Calculate the relative humidity from the vapour pressure and temperature.

    Args:
        p_v: vapour pressure, Pa
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        relative humidity, %
    """

    return p_v / svp.svp(equation=equation, status=status, t=t) * 100


def rh_to_pv(phi: float, t: float, equation: str, status: str) -> float:
    """Calculate the vapour pressure from the relative humidity and temperature

    Args:
        phi: relative humidity, %
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        vapour pressure, Pa
    """

    return svp.svp(equation=equation, status=status, t=t) * phi / 100


def rh_to_ah(phi: float, t: float, equation: str, status: str) -> float:
    """Calculate absolute humidity from relative humidity and temperature

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
    """Calculate relative humidity from absolute humidity and temperature

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


def t_dewp(x_s: float, equation: str, status: str) -> float:
    """Calculate the dew point temperature from the absolute humidity.

    Args:
        x_s: absolute humidity, kg/kgDA
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        dew point temperature, K
    """

    def f(t):
        return _get_saturated_absolute_humidity(t, equation, status) - x_s

    return scipy.optimize.newton(f, 273.15)


def _get_saturated_absolute_humidity(t: float, equation: str, status: str) -> float:
    """Calculate the saturated absolute humidity from the temperature

    Args:
        t: temperature, K
        equation: used equation for saturated vapour pressure
        status: ice or water

    Returns:
        saturated absolute humidity, kg/kgDA
    """

    p_sv = svp.svp(equation, status, t)

    return pv_to_ah(p_sv)


