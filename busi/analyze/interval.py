from enum import Enum
import numpy as np

"""
Time interval module
"""


class Interval(Enum):
    """the interval of the time

    The intervals described below are prepared.
        1 hour
        30 minutes
        15 minutes
    Notes:
        1時間を正数で分割する必要があるため、列挙型で管理することにした。
        上記に挙げたインターバール以外にも例えば以下のインターバールを将来的に実装することが考えられる。
        20分, 12分, 10分, 6分, 5分, 4分, 3分, 2分, 1分
    """

    H1 = '1h'
    M30 = '30m'
    M15 = '15m'

    def get_n_hour(self) -> int:
        """Calculate the number of steps between one hour.
        
        Returns:
            the number of steps between one hour

        Notes:
            1 hour: 1
            30 minutes: 2
            15 minutes: 4
        """

        return {
            Interval.H1: 1,
            Interval.M30: 2,
            Interval.M15: 4
        }[self]

    def get_time(self) -> float:
        """Get the interval time(hour).
        
        Returns:
            the interval time(hour), h
        """

        return {
            Interval.H1: 1.0,
            Interval.M30: 0.5,
            Interval.M15: 0.25
        }[self]

    def get_delta_t(self) -> int:
        """Get the interval time(seconds).
        
        Returns:
            the interval time(seconds), s
        """

        return {
            Interval.H1: 3600,
            Interval.M30: 1800,
            Interval.M15: 900
        }[self]

    def get_annual_number(self) -> int:
        """Get the number of the steps in a year.
        
        Returns:
            the number of the steps in a year
        Notes:
            If the result of the instanteneous value is treated,
            it may be needed that the value 1/1 0:00 and 12/31 24:00(=1/1 0:00 of the next year) should be prepared.
            In this case, One should be plused to the return value of this function.
            The number of this function is NOT include 12/31 24:00 (not one plus value).
        """

        return 8760 * self.get_n_hour()

    def get_pandas_freq(self) -> str:
        """Get the string argument of the frequency for pandas.
        
        Returns:
            freq argument of pandas
        """

        return {
            Interval.M15: '15min',
            Interval.M30: '30min',
            Interval.H1: 'H'
        }[self]

    def get_d_ns(self) -> np.ndarray:
        """Get the days of the year at step n. 1/1 is 1.

        Returns:
            the sequence of the days of the year at step n, d, [n]
        Notes:
            1/1 is 1.
            The return value is Ndarray of one dimension of the length (365 * 24 * n_hour).
            The value n_hour is the step of the hour.
                1h: 1
                30m: 2
                15m: 4
            Example of the output (in the case n_hour =1)
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,.....365,365,365,365
        """

        # the number of the step of one hour
        n_hour = self.get_n_hour()

        return np.repeat(np.arange(365) + 1, 24 * n_hour)

    def get_t_m_ns(self) -> np.ndarray:
        """Calculate the standard time at step n.

        Returns:
            the sequence of the standard time at step n, d, [n]
        """

        # the number of divided in one hour
        n_hour = self.get_n_hour()

        # the interval time, h
        int_interval = self.get_time()

        # 1h: 0, 1.0, .... , 23.0, 0, 1.0, ...23.0
        # 30m: 0, 0.5, 1.0, 1.5, .... , 23.5, 0, 0.5, ...23.5
        # 15m: 0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, .... , 23.75, 0, 0.25, ...23.75
        return np.tile(np.arange(24 * n_hour) * int_interval, 365)


