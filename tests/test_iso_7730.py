import unittest

from src.iso_7730 import get_pmv_ppd


class TestPMV_PPD(unittest.TestCase):

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


if __name__ == "__main__":
    unittest.main()
