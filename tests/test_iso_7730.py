import unittest

from src.iso_7730 import get_pmv_ppd, get_p_a


class TestPMV_PPD(unittest.TestCase):

    def test_get_p_a(self):

        expected = get_p_a(rh=50.0, t_a=26.0)
        actual = 1681.93  # Pa
        self.assertAlmostEqual(actual, expected, delta=0.01)

    def test_get_pmv_ppd(self):
        test_patterns = [
            ('01', 22.0, 22.0, 0.10, 60.0, 1.2, 0.5, -0.75, 17.0),
            ('02', 27.0, 27.0, 0.10, 60.0, 1.2, 0.5, 0.77, 17.0),
            ('03', 27.0, 27.0, 0.30, 60.0, 1.2, 0.5, 0.44, 9.0),
            ('04', 23.5, 25.5, 0.10, 60.0, 1.2, 0.5, -0.01, 5.0),
            ('05', 23.5, 25.5, 0.30, 60.0, 1.2, 0.5, -0.55, 11.0),
            ('06', 19.0, 19.0, 0.10, 40.0, 1.2, 1.0, -0.60, 13.0),
#            ('07', 23.5, 23.5, 0.10, 40.0, 1.2, 1.0, 0.50, 10.0),
            # TODO: don't matched.
            ('08', 23.5, 23.5, 0.30, 40.0, 1.2, 1.0, 0.12, 5.0),
            ('09', 23.0, 21.0, 0.10, 40.0, 1.2, 1.0, 0.05, 5.0),
            ('10', 23.0, 21.0, 0.30, 40.0, 1.2, 1.0, -0.16, 6.0),
            ('11', 22.0, 22.0, 0.10, 60.0, 1.6, 0.5, 0.05, 5.0),
            ('12', 27.0, 27.0, 0.10, 60.0, 1.6, 0.5, 1.17, 34.0),
            ('13', 27.0, 27.0, 0.30, 60.0, 1.6, 0.5, 0.95, 24.0)
        ]

        for run_num, t_a, t_r_bar, v_ar, rh, met_value, clo_value, expected_pmv, expected_ppd in test_patterns:
            with self.subTest(run_num=run_num):
                actual_pmv, actual_ppd = get_pmv_ppd(
                    met_value=met_value, p_eff=0.0, t_a=t_a, t_r_bar=t_r_bar, clo_value=clo_value, v_ar=v_ar, rh=rh)
                self.assertAlmostEqual(actual_pmv, expected_pmv, delta=0.01)
                self.assertAlmostEqual(actual_ppd, expected_ppd, delta=0.5)

    def test_get_pmv_ppd2(self):

        test_patterns = [
            ('01', 15.0, 15.0, 0.20, 50.0, 1.1, 0.8, -2.66, 96.00),
            ('02', 20.0, 15.0, 0.20, 50.0, 1.1, 0.8, -1.80, 67.00),
            ('03', 15.0, 20.0, 0.20, 50.0, 1.1, 0.8, -2.11, 82.00),
            ('04', 15.0, 15.0, 0.20, 90.0, 1.1, 0.8, -2.50, 93.00),
            ('05', 15.0, 15.0, 1.50, 50.0, 1.1, 0.8, -4.05, 100.00),
            ('06', 15.0, 15.0, 0.20, 50.0, 1.1, 1.1, -1.81, 68.00),
            ('07', 15.0, 15.0, 0.20, 50.0, 2.0, 0.8, -0.43, 9.14),
        ]

        for run_num, t_a, t_r_bar, v_ar, rh, met_value, clo_value, expected_pmv, expected_ppd in test_patterns:
            with self.subTest(run_num=run_num):
                actual_pmv, actual_ppd = get_pmv_ppd(
                    met_value=met_value, p_eff=0.0, t_a=t_a, t_r_bar=t_r_bar, clo_value=clo_value, v_ar=v_ar, rh=rh)
                self.assertAlmostEqual(actual_pmv, expected_pmv, delta=0.01)
                self.assertAlmostEqual(actual_ppd, expected_ppd, delta=0.5)


if __name__ == "__main__":
    unittest.main()
