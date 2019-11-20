import unittest

import src.vapor_pressure as vp


class TestVapourPressure(unittest.TestCase):

    def test_pv_to_ah(self):

        actual = vp.pv_to_ah(500.0)
        expected = 0.0030845524423506075
        self.assertAlmostEqual(actual, expected, delta=expected * 0.0001)

    def test_ah_to_pv(self):

        actual = vp.ah_to_pv(0.0030845524423506075)
        expected = 500.0
        self.assertAlmostEqual(actual, expected, delta=expected * 0.0001)

    def test_pv_to_rh(self):

        actual = vp.pv_to_rh(1000, 303.15, 'SONNTAG', 'water')
        expected = 23.54587078514155
        self.assertAlmostEqual(actual, expected, delta=expected * 0.0001)

    def test_rh_to_pv(self):

        actual = vp.rh_to_pv(23.54587078514155, 303.15,'SONNTAG','water')
        expected = 1000.0
        self.assertAlmostEqual(actual, expected, delta=expected * 0.0001)

    def test_rh_to_ah(self):

        actual = vp.rh_to_ah(50.0, 298.15, 'SONNTAG', 'water')
        expected = 0.009884095183878179
        self.assertAlmostEqual(actual, expected, delta=expected * 0.0001)

    def test_ah_to_rh(self):

        actual = vp.ah_to_rh(0.009884095183878179, 298.15, 'SONNTAG', 'water')
        expected = 50.0
        self.assertAlmostEqual(actual, expected, delta=expected * 0.0001)

    def test_get_dew_point_temp(self):

        actual = vp.get_dew_point_temp(0.014699216458052453, 'SONNTAG', 'water')
        expected = 293.15
        self.assertAlmostEqual(actual, expected, delta=expected * 0.0001)


if __name__ == "__main__":
    unittest.main()
