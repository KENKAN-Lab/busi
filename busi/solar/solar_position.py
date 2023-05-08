import math
from typing import Tuple, Optional


def calc(phi_loc: float, lambda_loc: float, tz: float, d: int, t_m: float) -> Tuple[float, float, float, float]:
    """Calculate the solar position, solar dicrination and the solar hour angle.

    Args:
        phi_loc: the latitude, rad
        lambda_loc: the longitude, rad
        tz: the time zone, the actual(clock) time for the location compared to UTC
        d: the days of the year
        t_m_ns: the standard hour

        solar declination, degrees
        the solar hour angle, degrees

    Returns:
        the solar altitude angle, the angle between the solar beam and the horizontal surface, rad
        the solar azimuth angle (angle from South, westwards positive, eastwards negative), rad
        the solar declination, rad
        the solar hour angle, rad
    """

    # 標準子午線(meridian), rad
    lambda_loc_mer = _get_lambda_loc_mer(tz=tz)

    # 1968年との年差
    n = _get_n()

    # 平均軌道上の近日点通過日（暦表時による1968年1月1日正午基準の日差）, d
    d_0 = _get_d_0(n=n)

    # 平均近点離角, rad
    m = _get_m(d=d, d_0=d_0)

    # 近日点と冬至点の角度, rad
    epsilon = _get_epsilon(m=m, n=n)

    # 真近点離角, rad
    v = _get_v(m=m)

    # 均時差, rad
    e_t = _get_e_t(m=m, epsilon=epsilon, v=v)

    # 赤緯, rad
    delta = _get_delta(epsilon=epsilon, v=v)

    # 時角, rad
    omega = _get_omega(t_m=t_m, lambda_loc=lambda_loc, lambda_loc_mer=lambda_loc_mer, e_t=e_t)

    # 太陽高度, rad
    h_sun = _get_h_sun(phi_loc=phi_loc, omega=omega, delta=delta)

    # 太陽の位置が天頂にないか（天頂にある = False, 天頂にない = True）
    is_not_zenith = _get_is_not_zenith(h_sun=h_sun)

    # 太陽の方位角の正弦（太陽が天頂に無い場合のみに定義される）
    sin_a_sun = _get_sin_a_sun(delta=delta, h_sun=h_sun, omega=omega, inzs=is_not_zenith)

    # 太陽の方位角の余弦（太陽が天頂に無い場合のみに定義される）
    cos_a_sun = _get_cos_a_sun(delta=delta, h_sun=h_sun, phi_loc=phi_loc, inzs=is_not_zenith)

    # 太陽の方位角（太陽が天頂に無い場合のみに定義される）, rad
    a_sun = _get_a_sun(cos_a_sun=cos_a_sun, sin_a_sun=sin_a_sun, inzs=is_not_zenith)

    return h_sun, a_sun, delta, omega


def _get_lambda_loc_mer(tz: float) -> float:
    """標準子午線を取得する。
    
    Args:
        the time zone, the actual(clock) time for the location compared to UTC

    Returns:
        標準子午線における経度, rad
    """

    return math.radians(15.0 * tz)


def _get_n() -> int:
    """Calculate the difference between the calculated year(1989) and 1968.

    Returns:
        the difference between the calculated year(1989) and 1968, year
    """

    # The solar position is calculated hear at the year 1989.
    y = 1989

    return y - 1968


def _get_d_0(n: int) -> float:
    """平均軌道上の近日点通過日を取得する。

    Args:
        n: 1968年との年差, year

    Returns:
        平均軌道上の近日点通過日（暦表時による1968年1月1日正午基準の日差）, d
    """

    return 3.71 + 0.2596 * n - int((n + 3.0) / 4.0)


def _get_m(d: int, d_0: float) -> float:
    """平均近点離角を計算する。
    
    Args:
        d: the day of the year (1/1 is 1), d
        d_0: 平均軌道上の近日点通過日（暦表時による1968年1月1日正午基準の日差）, d

    Returns:
        平均近点離角, rad
    """

    # 近点年（近日点基準の公転周期日数）
    d_ay = 365.2596

    # 平均近点離角, rad
    m_ns = 2 * math.pi * (d - d_0) / d_ay

    return m_ns


def _get_epsilon(m: float, n: int) -> float:
    """近日点と冬至点の角度を計算する。

    Args:
        m: 平均近点離角, rad
        n: 1968年との年差

    Returns:
        近日点と冬至点の角度, rad
    """

    return math.radians(12.3901 + 0.0172 * (n + m / (2 * math.pi)))


def _get_v(m: float) -> float:
    """真近点離角を計算する。

    Args:
        m: 平均近点離角, rad

    Returns:
        真近点離角, rad
    """

    return m + math.radians(1.914 * math.sin(m) + 0.02 * math.sin(2 * m))


def _get_e_t(m: float, epsilon: float, v: float) -> float:
    """

    Args:
        m: 平均近点離角, rad
        epsilon: 近日点と冬至点の角度, rad
        v: 真近点離角, rad
    Returns:
        均時差, rad
    """

    e_t = (m - v) - math.atan(0.043 * math.sin(2.0 * (v + epsilon)) / (1.0 - 0.043 * math.cos(2.0 * (v + epsilon))))

    return e_t


def _get_delta(epsilon: float, v: float) -> float:
    """赤緯を計算する。

    Args:
        epsilon: 近日点と冬至点の角度, rad
        v: 真近点離角, rad

    Returns:
        赤緯, rad

    Notes:
        赤緯は -π/2 ～ 0 π/2 の値をとる
    """

    # 北半球の冬至の日赤緯, rad
    delta_0 = math.radians(-23.4393)

    # 赤緯, rad
    delta = math.asin(math.cos(v + epsilon) * math.sin(delta_0))

    return delta


def _get_omega(t_m: float, lambda_loc: float, lambda_loc_mer: float, e_t: float) -> float:
    """時角を計算する。

    Args:
        t_m: the standard hour, h
        lambda_loc: the longitude, rad
        lambda_loc_mer: the longitude of the location of the standard hour, rad
        e_t: 均時差, rad

    Returns:
        時角, rad
    """

    return math.radians((t_m - 12.0) * 15.0) + (lambda_loc - lambda_loc_mer) + e_t


def _get_h_sun(phi_loc: float, omega: float, delta: float) -> float:
    """太陽高度を計算する。

    Args:
        phi_loc: the longitude, rad
        omega: 時角, rad
        delta: 赤緯, rad

    Returns:
        太陽高度, rad

    Notes:
        太陽高度はマイナスの値もとり得る。（太陽が沈んでいる場合）
    """

    h_sun = math.asin(math.sin(phi_loc) * math.sin(delta) + math.cos(phi_loc) * math.cos(delta) * math.cos(omega))

    return h_sun


def _get_is_not_zenith(h_sun: float) -> bool:
    """

    Args:
        h_sun_ns: 太陽高度, rad

    Returns:
        太陽の位置が天頂にないか（天頂にある = False, 天頂にない = True）
    """

    return h_sun != math.pi / 2


def _get_sin_a_sun(delta: float, h_sun: float, omega: float, inzs: bool) -> Optional[float]:
    """

    Args:
        delta: 赤緯, rad
        h_sun: 太陽高度, rad
        omega: 時角, rad
        inzs: 太陽位置が天頂にあるか否か（True=天頂にない, False=天頂にある）

    Returns:
        太陽の方位角の正弦（太陽が天頂に無い場合のみに定義される）
    """

    if inzs:
        return math.cos(delta) * math.sin(omega) / math.cos(h_sun)
    else:
        return None



def _get_cos_a_sun(delta: float, h_sun: float, phi_loc: float, inzs: bool) -> Optional[float]:
    """

    Args:
        delta: 赤緯, rad
        h_sun: 太陽高度, rad
        phi_loc: 緯度, rad
        inzs: 太陽位置が天頂にあるか否か（True=天頂にない, False=天頂にある）

    Returns:
        太陽の方位角の余弦（太陽が天頂に無い場合のみに定義される）
    """

    if inzs:
        return (math.sin(h_sun) * math.sin(phi_loc) - math.sin(delta)) / (math.cos(h_sun) * math.cos(phi_loc))
    else:
        return None


def _get_a_sun(cos_a_sun: float, sin_a_sun: float, inzs: bool) -> Optional[float]:
    """

    Args:
        cos_a_sun: 太陽の方位角の余弦（太陽が天頂に無い場合のみ計算する）
        sin_a_sun: 太陽の方位角の正弦（太陽が天頂に無い場合のみ計算する）
        inzs: 太陽位置が天頂にあるか否か（True=天頂にない, False=天頂にある）

    Returns:
        太陽の方位角（太陽が天頂に無い場合のみに定義される）, rad
    """

    if inzs:
        return math.atan2(sin_a_sun, cos_a_sun)
    else:
        return None
