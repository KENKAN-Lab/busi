import math
from collections import namedtuple
from typing import Tuple


def svp(equation: str, status: str, t: float) -> float:
    """Calculate the saturated vapor pressure.

    Args:
        equation: used equation name
        status: 'water' or 'ice'
        t: temperature, K

    Returns:
        saturated vapor pressure, Pa
    
    Notes:
        The equations prepared are below and the equations with asterisk need the status of water or ice.
        - SONNTAG *
        - WMO
        - WexlerHyland *
        - Tetens *
        - BriggsSacket *
        - Antoine *
        - GoffGratch *
    """

    if t < 0:
        raise ValueError('ERROR: Temperature can not be less than 0 K.')

    if equation == 'SONNTAG':
        return _saturated_vapor_pressure_SONNTAG(status, t)[0]
    elif equation == 'WMO':
        return _saturated_vapor_pressure_WMO(t)[0]
    elif equation == 'WexlerHyland':
        return _saturated_vapor_pressure_WH(status, t)[0]
    elif equation == 'Tetens':
        return _saturated_vapor_pressure_tetens(status, t)[0]
    elif equation == 'BriggsSacket':
        return _saturated_vapor_pressure_BS(status, t)[0]
    elif equation == 'Antoine':
        return _saturated_vapor_pressure_antoine(t)[0]
    elif equation == 'GoffGratch':
        return _saturated_vapor_pressure_GoffGratch(status, t)[0]
    else:
        raise ValueError('ERROR: false name of saturated vapor equation')


def dsvp_dt(equation: str, status: str, t: float) -> float:
    """Calculate the differential of the saturated vapor pressure.

    Args:
        equation: used equation name
        status: 'water' or 'ice'
        t: temperature, K

    Returns:
        differential of saturated vapor pressure, Pa/K

    Notes:
        The equations prepared are below and the equations with asterisk need the status of water or ice.
        - SONNTAG *
        - WMO
        - WexlerHyland *
        - Tetens *
        - BriggsSacket *
        - Antoine *
        - GoffGratch *
    """

    if t < 0:
        raise ValueError('ERROR: Temperature can not be less than 0 K.')

    if equation == 'SONNTAG':
        return _saturated_vapor_pressure_SONNTAG(status, t)[1]
    elif equation == 'WMO':
        return _saturated_vapor_pressure_WMO(t)[1]
    elif equation == 'WexlerHyland':
        return _saturated_vapor_pressure_WH(status, t)[1]
    elif equation == 'Tetens':
        return _saturated_vapor_pressure_tetens(status, t)[1]
    elif equation == 'BriggsSacket':
        return _saturated_vapor_pressure_BS(status, t)[1]
    elif equation == 'Antoine':
        return _saturated_vapor_pressure_antoine(t)[1]
    elif equation == 'GoffGratch':
        return _saturated_vapor_pressure_GoffGratch(status, t)[1]
    else:
        raise ValueError('ERROR: false name of saturated vapor equation')


def _saturated_vapor_pressure_SONNTAG(status: str, t: float) -> Tuple[float, float]:
    """Calculate the saturated vapor pressure and its differential.

    Args:
        status: 'water' or 'ice'
        t: temperature, K

    Returns:
        2 parameters:
            (1) saturated vapor pressure, Pa
            (2) differential of saturated vapor pressure, Pa/K
    """

    Coeff = namedtuple('Coeff', ('a1', 'a2', 'a3', 'a4', 'a5'))

    c = {
        'water': Coeff(-6096.9385, 21.2409642, -0.02711193, 0.00001673952, 2.433502),
        'ice': Coeff(-6024.5282, 29.32707, 0.010613863, -0.000013198825, -0.49382577)
    }[status]

    k = c.a1 / t + c.a2 + c.a3 * t + c.a4 * t ** 2 + c.a5 * math.log(t)

    pvs = math.exp(k)

    dpvs_dt = pvs * (- c.a1 / (t ** 2) + c.a3 + 2 * c.a4 * t + c.a5 / t)

    return pvs, dpvs_dt


def _saturated_vapor_pressure_WMO(t: float) -> Tuple[float, float]:
    """Calculate the saturated vapor pressure and its differential.

    Args:
        t: temperature, K

    Returns:
        2 parameters:
            (1) saturated vapor pressure, Pa
            (2) differential of saturated vapor pressure, Pa/K
    """

    ew = 2.78614 + 10.79574 * (1.0 - 273.16 / t) \
        - 5.028 * math.log10(t / 273.16) \
        + 1.50475 * 10**(-4) * (1.0 - 10.0 ** (-8.2969 * (t / 273.16 - 1.0))) \
        + 0.42873 * 10**(-3) * (10 ** (4.76955 * (1.0 - 273.16 / t)) - 1.0)

    dew_dt = 10.79574 * 273.16 / (t ** 2) \
        - 5.028 / math.log(10) / t \
        + 1.50475 * 10 ** (-4) * 8.2969 / 273.16 * math.log(10.0) * 10.0 ** (-8.2969 * (t / 273.16 - 1.0)) \
        + 0.42873 * 10 ** (-3) * 4.76955 * 273.16 / (t ** 2) * math.log(10.0) * 10 ** (4.76955 * (1.0 - 273.16 / t))

    return 10**ew, 10**ew * dew_dt * math.log(10.0)


def _saturated_vapor_pressure_WH(status: str, t: float) -> Tuple[float, float]:
    """Calculate the saturated vapor pressure and its differential.

    Args:
        status: 'water' or 'ice'
        t: temperature, K

    Returns:
        2 parameters:
            (1) saturated vapor pressure, Pa
            (2) differential of saturated vapor pressure, Pa/K
    """

    Coeff = namedtuple('Coeff', ('a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'a4', 'b4', 'a5', 'b5', 'a6', 'b6', 'a7', 'b7'))

    c = {
        'water': Coeff(-0.58002206, 4, 0.13914993, 1, -0.48640239, -1, 0.41764768, -4, -0.14452093, -7,  0.0, 0, 0.65459673, 1),
        'ice': Coeff(-0.56745359, 4, 0.63925247, 1, -0.96778430, -2, 0.62215701, -6,  0.20747825, -8, -0.94840240, -12, 0.41635019, 1)
    }[status]

    k = c.a1 * 10 ** c.b1 / t + c.a2 * 10 ** c.b2 + c.a3 * 10 ** c.b3 * t + c.a4 * 10 ** c.b4 * t ** 2 + c.a5 * 10 ** c.b5 * t ** 3 + c.a6 * 10 ** c.b6 * t ** 4 + c.a7 * 10 ** c.b7 * math.log(t)

    pvs = math.exp(k)

    dpvs_dt = pvs * (-c.a1 * 10 ** c.b1 / t ** 2 + c.a3 * 10 ** c.b3 + 2 * c.a4 * 10 ** c.b4 * t + 3 * c.a5 * 10 ** c.b5 * t ** 2 + 4 * c.a6 * 10 ** c.b6 * t ** 3 + c.a7 * 10 ** c.b7 / t)

    return pvs, dpvs_dt


def _saturated_vapor_pressure_tetens(status: str, t: float) -> Tuple[float, float]:
    """Calculate the saturated vapor pressure and its differential.

    Args:
        status: 'water' or 'ice'
        t: temperature, K

    Returns:
        2 parameters:
            (1) saturated vapor pressure, Pa
            (2) differential of saturated vapor pressure, Pa/K
    """

    Coeff = namedtuple('Coeff', ('a', 'b'))
    c = {
        'water': Coeff(17.2693882, 35.86),
        'ice': Coeff(21.8745584, 7.66)
    }[status]

    pvs = 6.1078 * 10**2 * math.exp(c.a * (t - 273.16) / (t - c.b))

    dpvs_dt = pvs * (c.a / (t - c.b) - c.a * (t - 273.16) / ((t - c.b) ** 2))

    return pvs, dpvs_dt


def _saturated_vapor_pressure_BS(status: str, t: float) -> Tuple[float, float]:
    """Calculate the saturated vapor pressure and its differential.

    Args:
        status: 'water' or 'ice'
        t: temperature, K

    Returns:
        2 parameters:
            (1) saturated vapor pressure, Pa
            (2) differential of saturated vapor pressure, Pa/K
    """

    Coeff = namedtuple('Coeff', ('a1', 'a2', 'a3', 'a4', 'a5'))

    c = {
        'water': Coeff(-2313.0338, -164.03307, 38.053682, -0.13844344, 0.000074465367),
        'ice': Coeff(-5631.1206, -8.363602, 8.2312, -0.03861449, 0.0000277494)
    }[status]

    k = c.a1 / t + c.a2 + c.a3 * math.log(t) + c.a4 * t + c.a5 * t ** 2 - math.log(10.0)

    pvs = math.exp(k)

    dpvs_dt = pvs * (-c.a1 / t ** 2 + c.a3 / t + c.a4 + 2 * c.a5 * t)

    return pvs, dpvs_dt


def _saturated_vapor_pressure_antoine(t: float) -> Tuple[float, float]:
    """Calculate the saturated vapor pressure and its differential.

    Args:
        t: temperature, K

    Returns:
        2 parameters:
            (1) saturated vapor pressure, Pa
            (2) differential of saturated vapor pressure, Pa/K
    """

    a = 8.02754
    b = 1705.616
    c = 231.405 - 273.15

    pvs = 10**(a - (b / (t + c)))

    dpvs_dt = pvs * math.log(10) * (b / (t + c) ** 2)

    return pvs * 101325 / 760, dpvs_dt * 101325 / 760  # 760mmHg = 101325 Pa


def _saturated_vapor_pressure_GoffGratch(status: str, t: float) -> Tuple[float, float]:
    """calculate the saturated vapor pressure and its differential

    Args:
        status: 'water' or 'ice'
        t: temperature, K

    Returns:
        2 parameters:
            (1) saturated vapor pressure, Pa
            (2) differential of saturated vapor pressure, Pa/K
    """

    a_1 = -7.90298
    a_2 = 5.02808
    a_3 = -1.3816 * 10**(-7)
    a_4 = 11.344
    a_5 = 8.1328 * 10**(-3)
    a_6 = -3.49149
    t_st = 373.15
    e_st = 1013.25
    k_w = a_1 * (t_st / t - 1) + a_2 * math.log10(t_st / t) + a_3 * (10 ** (a_4 * (1 - t / t_st)) - 1) + a_5 * (10 ** (a_6 * (t_st / t - 1)) - 1) + math.log10(e_st)
    pvs_w = 10**k_w
    dpvs_dt_w = pvs_w * math.log(10) * (- a_1 * t_st / t ** 2 - a_2 / t / math.log(10) - a_3 * math.log(10) * 10 ** (a_4 * (1 - t / t_st)) * a_4 / t_st - a_5 * math.log(10) * 10 ** (a_6 * (t_st / t - 1)) * a_6 * t_st / t ** 2)

    b_1 = -9.09718
    b_2 = -3.56654
    b_3 = 0.876793
    T_0 = 273.16
    e_i0 = 6.1173
    k_i = b_1 * (T_0 / t - 1) + b_2 * math.log10(T_0 / t) + b_3 * (1 - t / T_0) + math.log10(e_i0)
    pvs_i = 10**k_i
    dpvs_dt_i = pvs_i * math.log(10) * (-b_1 * T_0 / t ** 2 - b_2 / math.log(10) / t - b_3 / T_0)

    return {
        'water': (pvs_w * 100, dpvs_dt_w * 100),
        'ice': (pvs_i * 100, dpvs_dt_i * 100)
    }[status]
