from enum import Enum, auto
import dataclasses
from typing import List, Dict
import numpy as np
from scipy import optimize


class VentilationType(Enum):

    TYPE1 = auto()
    TYPE2 = auto()
    TYPE3 = auto()


@dataclasses.dataclass
class Gap:

    # height of gap, m
    height: float

    # ratio of gap, 0.0 ~ 1.0
    ratio: float


class Port:

    def __init__(self, height: float, is_specified: bool, oa: float = 0.0, ratio: float = 0.0):
        """

        Args:
            height: height of the port, m
            is_specified: the equivalent opening area is specified ?
            oa: equivalent opening area, m2
            ratio: ratio of the port area (0.0 ~ 1.0)
        """
        self._height = height
        self._is_specified = is_specified
        self._oa = oa
        self._ratio = ratio

    def get_equivalent_opening_area(self, v_vent: float, vent_type: VentilationType):

        if vent_type == VentilationType.TYPE1:
            return 0.0
        else:
            if self._is_specified:
                return self._oa
            else:
                return v_vent * 0.7 / 10000 * self._ratio

    @property
    def height(self):
        return self._height


class SimpleVentilation:

    def __init__(
            self,
            c_value: float,
            a_total: float,
            vent_type: VentilationType,
            gaps: List[Gap],
            ports: List[Port],
            h_room: float = 2.4
    ):
        """

        Args:
            c_value: C value, cm2/m2
            a_total: total floor area, m2
            vent_type: VentilationType
            gaps: list of gap
            ports: list of ports
            h_room: room height, m
        """

        self._c_value = c_value
        self._a_total = a_total
        self._vent_type = vent_type
        self._gaps = gaps
        self._ports = ports
        self._h_room = h_room

        # room volume, m3
        self._v_room = self._a_total * self._h_room

        # ventilation amount, m3/h
        self._v_vent = self._v_room * 0.5

        # height of opening, m [n]
        self._hs = np.concatenate([
            np.array([g.height for g in self._gaps]),
            np.array([p.height for p in self._ports])
        ])

        # equivalent opening area, m2 [n]
        oas_gap = np.array([self._c_value * self._a_total * g.ratio / 10000 for g in self._gaps])
        oas_port = np.array(
            [p.get_equivalent_opening_area(v_vent=self._v_vent, vent_type=self._vent_type) for p in self._ports]
        )
        self._oas = np.concatenate([oas_gap, oas_port])

    def get_vent(self, theta_in: float, theta_out: float, vent_ratio: float) -> Dict:
        """

        Args:
            theta_in: inside temperature, degree C
            theta_out: outside temperature, degree C
            vent_ratio: the ratio of actual ventilation amount for the 0.5 frequency amount, (0.0 ~ 1.0)

        Returns:

        """

        # density, kg/m3
        rho_in = self.get_density(theta_in)
        rho_out = self.get_density(theta_out)

        def get_q_total(p0: float):
            """
            Args:
                p0: pressure difference (positive: negative pressure), Pa
            Returns:
                total balance of ventilation (outside to inside), m3/h
            """

            return np.sum(self.get_qs(p0=p0, rho_out=rho_out, rho_in=rho_in)) + {
                VentilationType.TYPE1: 0.0,
                VentilationType.TYPE2: self._v_vent * vent_ratio,
                VentilationType.TYPE3: -self._v_vent * vent_ratio
            }[self._vent_type]

        # level 0 pressure, Pa
        p0 = float(optimize.fsolve(func=get_q_total, x0=np.array([0.0]))[0])
        print(p0)

        # amount of ventilation (outside to inside), m3/h, [n]
        qs = self.get_qs(p0=p0, rho_out=rho_out, rho_in=rho_in)

        ps = self.get_ps(p0=p0, rho_out=rho_out, rho_in=rho_in)

        h_n = self.get_natural_band_height(p0=p0, rho_out=rho_out, rho_in=rho_in)

        gap_vent = {
            VentilationType.TYPE1: np.sum(qs[qs > 0.0]),
            VentilationType.TYPE2: np.sum(qs[qs > 0.0]),
            VentilationType.TYPE3: np.sum(qs[qs < 0.0]),
        }[self._vent_type]

        gap_vent_time = gap_vent / self._v_room

        return {
            'p0': p0,
            'qs': qs,
            'ps': ps,
            'h_n': h_n,
            'gap_vent': gap_vent,
            'gap_vent_time': gap_vent_time
        }

    @staticmethod
    def get_natural_band_height(p0: float, rho_out: float, rho_in: float) -> float:
        """
        get natural band height
        Args:
            p0: pressure of the 0 height point, Pa
            rho_out: outside air density, kg/m3
            rho_in: inside air density, kg/m3

        Returns:
            natural band height, m
        """

        return p0 / (9.8 * (rho_out - rho_in))

    def get_qs(self, p0: float, rho_out: float, rho_in: float) -> np.ndarray:
        """
        get pressure difference between inside and outside
        Args:
            p0: pressure of the 0 height point, Pa
            rho_out: outside air density, kg/m3
            rho_in: inside air density, kg/m3

        Returns:
            amount of ventilation (outside to inside), m3/h, [n]
        """

        ps = self.get_ps(p0=p0, rho_out=rho_out, rho_in=rho_in)

        return np.sign(ps) * self._oas * (2.0 / np.where(ps > 0.0, rho_out, rho_in) * np.abs(ps)) ** 0.5 * 3600

    def get_ps(self, p0: float, rho_out: float, rho_in: float) -> np.ndarray:
        """
        get pressure difference between inside and outside
        Args:
            p0: pressure of the 0 height point, Pa
            rho_out: outside air density, kg/m3
            rho_in: inside air density, kg/m3

        Returns:
            pressure difference between inside and outside, Pa, [n]
        """

        return p0 - 9.8 * self._hs * (rho_out - rho_in)

    @staticmethod
    def get_density(theta: float):
        """
        get air density
        Args:
            theta: temperature, degree C

        Returns:
            density, kg/m3
        """

        return 353.25 / (theta + 273.15)


def get_simple_calc(story: int, theta_in: float, theta_out: float, c_value: float, a_total: float, vent_type: VentilationType):

    gaps = {
        1: [
            Gap(height=0.0, ratio=0.5),
            Gap(height=2.4, ratio=0.5)
        ],
        2: [
            Gap(height=0.0, ratio=0.25),
            Gap(height=2.4, ratio=0.25),
            Gap(height=2.9, ratio=0.25),
            Gap(height=5.3, ratio=0.25)
        ]
    }[story]

    ports = {
        1: [
            Port(height=1.6, is_specified=False, ratio=1.0)
        ],
        2: [
            Port(height=1.6, is_specified=False, ratio=0.5),
            Port(height=4.5, is_specified=False, ratio=0.5)
        ]
    }[story]

    v = SimpleVentilation(c_value=c_value, a_total=a_total, gaps=gaps, ports=ports, vent_type=vent_type)

    return v.get_vent(theta_in=theta_in, theta_out=theta_out, vent_ratio=1.0)


if __name__ == '__main__':

    result = get_simple_calc(
        story=2,
        theta_in=20.0,
        theta_out=0.0,
        c_value=2.0,
        a_total=120.0,
        vent_type=VentilationType.TYPE1
    )

    print(result)
