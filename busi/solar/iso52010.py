# This module is based on ISO 52010 part 1
# "Energy performance of buildings -External climatic conditions- Part 1 Conversion of climatic data for energy calculations"

import numpy as np
import math
from typing import Tuple


def calc(n_day: int, tz: float, lambda_w: float, n_hour: float, phi_w: float, half_hour_shift: bool = True) -> Tuple[float, float, float, float]:
    """Calculate the solar position, the solar diclination and the solar hour angle.

    Args:
        n_day: the day of the year, from 1 to 365 or 366 (leap year)
        tz: the time zone, the actual(clock) time for the location compared to UTC
        lambda_w: the longitude of the weather station, degrees
        n_hour: the actual(clock) time for the location, the hour of the day
        phi_w: the latitude of the weather station, degrees
        half_hour_sift: does the 30 min shift? (Please see the Note of the function "_get_omega".)

    Returns:
        the solar altitude angle, the angle between the solar beam and the horizontal surface, degrees
        the solar azimuth angle (angle from South, eastwards positive, westwards negative), degrees
        solar declination, degrees
        the solar hour angle, degrees
    """

    delta = _get_delta(n_day)
    t_eq = _get_t_eq(n_day)
    t_shift = _get_t_shift(tz, lambda_w)
    t_sol = _get_t_sol(n_hour, t_eq, t_shift)
    omega = _get_omega(t_sol, half_hour_shift)
    alpha_sol = _get_alpha_sol(delta, phi_w, omega)
    phi_sol = _get_phi_sol(delta, omega, alpha_sol, phi_w)

    return alpha_sol, phi_sol, delta, omega


def _get_delta(n_day: int) -> float:
    """Calculate the solar diclination.    
    
    Args:
        n_day: the day of the year, from 1 to 365 or 366 (leap year)
    
    Returns:
        solar declination, degrees
    """

    # the earth orbit deviation, degrees
    r_dc = 360 / 365 * n_day

    delta = 0.33281 \
          - 22.984 * math.cos(math.radians(r_dc)) \
          - 0.3499 * math.cos(math.radians(2 * r_dc)) \
          - 0.1398 * math.cos(math.radians(3 * r_dc)) \
          + 3.7872 * math.sin(math.radians(r_dc)) \
          + 0.03205 * math.sin(math.radians(2 * r_dc)) \
          + 0.07187 * math.sin(math.radians(3 * r_dc))

    return delta


def _get_t_eq(n_day: int) -> float:
    """Calculate the equation of time.

    Args:
        n_day: the day of the year, from 1 to 365 or 366 (leap year)

    Returns:
        the equation of time, min
    """

    if n_day < 21:
        return 2.6 + 0.44 * n_day
    elif n_day < 136:
        return 5.2 + 0.0 * math.cos((n_day - 43) * 0.0357 * 180 / math.pi)
    elif n_day < 241:
        return 1.4 - 5.0 * math.cos((n_day - 135) * 0.0449 * 180 / math.pi)
    elif n_day < 336:
        return -6.3 - 10.0 * math.cos((n_day - 306) * 0.036 * 180 / math.pi)
    else:
        return 0.45 * (n_day - 359)


def _get_t_shift(tz: float, lambda_w: float) -> float:
    """Calculate the time shift.

    Args:
        tz: the time zone, the actual(clock) time for the location compared to UTC
        lambda_w: the longitude of the weather station, degrees

    Returns:
        the time sift, h
    
    Note:
        Daylight saving time is sdisregarded in t_shift whidh is time independent.
    """

    return tz - lambda_w / 15


def _get_t_sol(n_hour: float, t_eq: float, t_shift: float) -> float:
    """Calculate the solar time.

    Args:
        n_hour: the actual(clock) time for the location, the hour of the day
        t_eq: the equation of time, min
        t_shift: the time shift, h

    Returns:
        the sloar time, h
    """

    return n_hour - t_eq / 60 - t_shift


def _get_omega(t_sol: float, half_hour_shift: bool=True) -> float:
    """Calculate the solar hour angle.

    Args:
        t_sol: the solar time, h
        half_hour_sift: does 30 min sift ? (Please see Note. The default value of ISO is true.)

    Returns:
        the solar hour angle, degrees
    
    Note:
        Explanation of "12.5":
        The hour numbers are actually hour sections:
        the first hour section of a day runs from 0h to 1h.
        So, the average position of the sun for the solar radiation measured during (solar) hour section N is at (solar) time = (N - 0.5) h of the (solar) day.
    """

    if half_hour_shift:
        omega = 180 / 12 * (12.5 - t_sol)
    else:
        omega = 180 / 12 * (12 - t_sol)

    if omega > 180.0:
        omega = omega - 360.0
    
    if omega < -180:
        omega = omega + 360.0
    
    return omega


def _get_alpha_sol(delta: float, phi_w: float, omega: float) -> float:
    """Calculate the solar altitude angle.

    The solar altitude angle is the angle between the solar beam and the horizontal surface,
    determined in the middle of the current hour as a function of the solar hour angle, the solar declination and the latitude.

    Args:
        delta: the solar declination, degrees
        phi_w: the latitude of the weather station, degrees
        omega: the solar hour angle of the weather station, degrees

    Returns:
        the solar altitude angle, the angle between the solar beam and the horizontal surface, degrees
    """

    alpha_sol = math.degrees(math.asin(
        math.sin(math.radians(delta)) * math.sin(math.radians(phi_w))
        + math.cos(math.radians(delta)) * math.cos(math.radians(phi_w)) * math.cos(math.radians(omega))
    ))

    if alpha_sol < 0.0001:
        alpha_sol = 0.0
    
    return alpha_sol


def _get_theta_z(alpha_sol: float) -> float:
    """Calculate the solar zenith angle, the angle between the solar beam and the zenith.

    Args:
        alpha_sol: the sloar altitude angle, the angle between the solar beam and the horizontal surface, degrees

    Returns:
        the solar zenith angle, the angle between the solar beam and the zenith, degrees
    """

    return 90.0 - alpha_sol


def _get_phi_sol(delta: float, omega: float, alpha_sol: float, phi_w: float) -> float:
    """Calculate the solar azimuth.

    Args:
        delta: the solar declination, degrees
        omega: the solar hour angle, degrees
        alpha_sol: the solar altitude angle, degrees
        phi_w: the latitude of the weather station, degrees
    
    Returns:
        the solar azimuth angle (angle from South, eastwards positive, westwards negative), degrees
    """

    sin_phi_sol_aux1 = (
        math.cos(math.radians(delta)) * math.sin(math.radians(180 - omega))
        / math.cos(math.asin(math.sin(math.radians(alpha_sol))))
    )

    cos_phi_sol_aux1 = (
        (
            math.cos(math.radians(phi_w)) * math.sin(math.radians(delta))
            + math.sin(math.radians(phi_w)) * math.cos(math.radians(delta)) * math.cos(math.radians(180 - omega))
        ) / math.cos(math.asin(math.sin(math.radians(alpha_sol))))
    )

    # The equation (15) in ISO 52010 part1 2017 is mistaken.
    phi_sol_aux2 = math.degrees(math.asin(
        math.cos(math.radians(delta)) * math.sin(math.radians(180 - omega))
        / math.cos(math.asin(math.sin(math.radians(alpha_sol))))
    ))
    
    if (sin_phi_sol_aux1 >= 0) and (cos_phi_sol_aux1 > 0):
        return 180 - phi_sol_aux2
    elif cos_phi_sol_aux1 < 0:
        return phi_sol_aux2
    else:
        return -(180 + phi_sol_aux2)


def _get_theta_sol_ic(delta: float, phi_w: float, beta_ic: float, gamma_ic: float, omega: float) -> float:
    """Calculate the solar angle of incidence on the inclined surface.

    Args:
        delta: the solar declination, degrees
        phi_w: the latitude of the weather station, degrees
        beta_ic: the tilt angle of the inclined surface, degrees
        gamma_ic: the orientation of the inclined surface ( angle from South, eastwards positive, westwards negative), degrees
        omega: the solar hour angle of the weather station, degrees

    Returns:
        the solar angle of incidence onthe inclined surface, degrees
    """

    return math.degrees(math.arccos(
          math.sin(math.degrees(delta)) * math.sin(math.degrees(phi_w)) * math.cos(math.radians(beta_ic))
        - math.sin(math.radians(delta)) * math.cos(math.radians(phi_w)) * math.sin(math.radians(beta_ic)) * math.cos(math.radians(gamma_ic))
        + math.cos(math.radians(delta)) * math.cos(math.radians(phi_w)) * math.cos(math.radians(beta_ic)) * math.cos(math.radians(omega))
        + math.cos(math.radians(delta)) * math.sin(math.radians(phi_w)) * math.sin(math.radians(beta_ic)) * math.cos(math.radians(gamma_ic)) * math.cos(math.radians(omega))
        + math.cos(math.radians(delta)) * math.sin(math.radians(beta_ic)) * math.sin(math.radians(gamma_ic)) * math.sin(math.radians(omega))
    ))


def _get_gamma_sol_ic(omega: float, gamma_ic: float) -> float:
    """Calculate the azimuth angle between sun and the inclined surface.

    Args:
        omega: the solar hour angle of the weather station, degrees
        gamma_ic: the orientation of the inclined surface, degrees

    Returns:
        the azimuth angle between sun and the inclined surface, degrees
    """

    if (omega - gamma_ic) < 180:
        return 360 + omega - gamma_ic
    else:
        return omega - gamma_ic


def _get_beta_sol_ic(beta_ic: float, theta_z: float) -> float:
    """Calculate the tilt angle between sun and inclined surface.

    Args:
        beta_ic: the tilt angle of the inclined surface, degrees
        theta_z: the solar zenith angle, the angle between the solar beam and the zenith, degrees

    Returns:
        the tilt angle between sun and the inclined surface, degrees
    """

    if (beta_ic - theta_z) > 180.0:
        return -360 + beta_ic - theta_z
    elif (beta_ic - theta_z) < -180.0:
        return 360 + beta_ic - theta_z
    else:
        return beta_ic - theta_z


def _get_m(alpha_sol: float) -> float:
    """Calculate the air mass.

    Args:
        alpha_sol: the solar altitude angle, degrees

    Returns:
        the dimensionless air mass
    """

    if alpha_sol >= 10.0:
        return 1 / math.sin(math.radians(alpha_sol))
    else:
        return 1 / (math.sin(math.radians(alpha_sol)) + 0.15 * (alpha_sol + 3.885)**-1.253)


