import math
from scipy.optimize import newton


def get_pmv_ppd(
        met_value: float, p_eff: float, t_a: float, t_r_bar: float, clo_value: float, v_ar: float, rh: float
) -> (float, float):
    """calculate PMV & PPD

    Args:
        t_a: the air temperature, degree C
        t_r_bar: the mean radiant temperature, degree C
        rh: the relative humidity, %
        v_ar: the relative air velocity, m/s
        met_value: the metabolic rate, met
        p_eff: the effective mechanical power, met
        clo_value: the clothing insulation, clo

    Returns:
        tuple:
            PMV
            PPD

    Notes:
        Reference: ISO 7730, 1994, 2005
    """

    # the water vapour partial pressure, Pa
    p_a = get_p_a(rh, t_a)

    # the clothing insulation, m2K/W
    i_cl = convert_clo_to_m2kw(clo_value)

    # the metabolic rate, W/m2
    m = convert_met_to_wm2(met_value)

    # the effective mechanical power, W/m2
    w = convert_met_to_wm2(p_eff)

    # the internal heat production in the human body, W/m2
    mw = m - w

    # the clothing surface area factor
    f_cl = get_f_cl(i_cl)

    t_cl = newton(lambda t: get_t_cl(f_cl, i_cl, mw, t_a, t, v_ar, t_r_bar) - t, 0.001)

    h_c = max(12.1 * math.sqrt(v_ar), 2.38 * abs(t_cl - t_a) ** 0.25)

    pmv = get_pmv(f_cl, h_c, m, p_a, t_a, t_cl, t_r_bar, w)

    return pmv


def get_t_cl(f_cl, i_cl, mw, t_a, t_cl, v_ar, t_r_bar):

    h_c = max(12.1 * math.sqrt(v_ar), 2.38 * abs(t_cl - t_a) ** 0.25)

    return 35.7 - 0.028 * mw - i_cl * (
                3.96 * 10 ** (-8) * f_cl * ((t_cl + 273.0) ** 4.0 - (t_r_bar + 273.0) ** 4.0)
                + f_cl * h_c * (t_cl - t_a)
    )


def get_pmv(f_cl, h_c, m, p_a, t_a, t_cl, t_r_bar, w):
    """

    Args:
        f_cl: the clothing insulation, m2K/W
        h_c: the convective heat transfer coefficient, W/m2K
        m: the metabolic rate, W/m2
        p_a: the water vapour partial pressure, Pa
        t_a: the air temperature, degree C
        t_cl: the clothing surface temperature, degree C
        t_r_bar: the mean radiant temperature, degree C
        w: the effective mechanical power, W/m2

    Returns:
        PMV

    Notes:
        equation (1)
    """

    pmv = (0.303 * math.exp(-0.036 * m) + 0.028) * (
            (m - w)
            - 3.05 * 10 ** (-3) * (5733.0 - 6.99 * (m - w) - p_a)
            - max(0.42 * ((m - w) - 58.15), 0.0)
            - 1.7 * 10 ** (-5) * m * (5867.0 - p_a)
            - 0.0014 * m * (34.0 - t_a)
            - 3.96 * 10 ** (-8) * f_cl * ((t_cl + 273) ** 4.0 - (t_r_bar + 273.0) ** 4.0)
            - f_cl * h_c * (t_cl - t_a))

    ppd = get_ppd(pmv=pmv)

    return pmv, ppd


def get_p_a(rh: float, t_a: float) -> float:
    """

    Args:
        rh: relative humidity, %
        t_a: the air temperature, degree C

    Returns:
        the water vapour partial pressure, Pa

    """

    # TODO 飽和水蒸気圧の計算方法は省エネ基準で採用している方法（WMO?）に揃えた方が良い。
    return rh / 100. * FNPS(t_a) * 1000.0


def convert_clo_to_m2kw(clo):
    """convert the unit of clo to m2K/W

    Args:
        clo: value, clo

    Returns:
        value, m2K/W
    """
    return clo * 0.155


def convert_met_to_wm2(met):
    """convert the unit of met to W/m2

    Args:
        met: value, met

    Returns:
        value, W/m2
    """

    return met * 58.15


def get_f_cl(i_cl: float) -> float:
    """calculate clothing surface area factor

    Args:
        i_cl: the clothing insulation, m2K/W

    Returns:
        the clothing surface area factor

    Notes:
        equation (4)
    """

    if i_cl <= 0.078:
        return 1.00 + 1.290 * i_cl
    else:
        return 1.05 + 0.645 * i_cl


def get_ppd(pmv: float) -> float:
    """calculate PPD

    Args:
        pmv: PMV

    Returns:
        PPD

    Notes:
        PPD is Predicted Percentage of Dissatisfied.
    """

    return 100.0 - 95.0 * math.exp(-0.03353 * pmv ** 4.0 - 0.2179 * pmv ** 2.0)



# 飽和水蒸気圧[kPa]の計算（ASHRAE Standard 55-2013）
def FNPS(T: float) -> float:
    return math.exp(16.6536 - 4030.183 / (T + 235.0))

