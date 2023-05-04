import unittest

from busi.window import multiple_reflection as mr


class GlassTest(unittest.TestCase):

    def test_something(self):

        ss = [
            mr.SolarSpecSingleLayer(tau_f=0.859, tau_b=0.859, rho_f=0.077, rho_b=0.077),
            mr.SolarSpecSingleLayer(tau_f=0.859, tau_b=0.859, rho_f=0.077, rho_b=0.077),
            mr.SolarSpecSingleLayer(tau_f=0.859, tau_b=0.859, rho_f=0.077, rho_b=0.077)
        ]

        g = mr.Glass(ss=ss)

        rho = g.get_total_solar_spec().rho_f
        tau = g.get_total_solar_spec().tau_f

        alpha = g.get_abs_multi_layer()

        # グレージング複合体全体の正面側の反射率・グレージング複合体全体の透過率・層jの日射吸収率の合計
        total = rho + tau + sum(alpha)

        self.assertAlmostEqual(total, 1.0)
        self.assertAlmostEqual(rho, 0.1770242008152192)
        self.assertAlmostEqual(tau, 0.6442755896857679)
        self.assertAlmostEqual(alpha[0], 0.0714523269524727)
        self.assertAlmostEqual(alpha[1], 0.05924597598089509)
        self.assertAlmostEqual(alpha[2], 0.04800190656564512)


if __name__ == '__main__':
    unittest.main()
