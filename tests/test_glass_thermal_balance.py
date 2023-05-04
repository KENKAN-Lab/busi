import unittest
from busi.window import glass_thermal_balance as t
import numpy as np


class GlassesTest(unittest.TestCase):

    def test_glass_unit(self):

        g1 = t.GlassUnit(d=0.003)
        self.assertEqual(g1.r, 0.003)

        g2 = t.GlassUnit(d=0.003, lmd=0.5)
        self.assertEqual(g2.r, 0.006)

    def test_surface(self):

        s1 = t.NonLowESurface()
        self.assertAlmostEqual(s1.eps, 0.837)

        s2 = t.LowESurface(epsilon_n=0.04)
        self.assertAlmostEqual(s2.eps, 0.048)

        s3 = t.LowESurface(epsilon_n=0.075)
        self.assertAlmostEqual(s3.eps, 0.087)

        s4 = t.LowESurface(epsilon_n=0.15)
        self.assertAlmostEqual(s4.eps, 0.168)

        s5 = t.LowESurface(epsilon_n=0.25)
        self.assertAlmostEqual(s5.eps, 0.27)

        s6 = t.LowESurface(epsilon_n=0.35)
        self.assertAlmostEqual(s6.eps, 0.36575)

        s7 = t.LowESurface(epsilon_n=0.45)
        self.assertAlmostEqual(s7.eps, 0.45675)

        s8 = t.LowESurface(epsilon_n=0.55)
        self.assertAlmostEqual(s8.eps, 0.5445)

        s9 = t.LowESurface(epsilon_n=0.65)
        self.assertAlmostEqual(s9.eps, 0.6305)

        s10 = t.LowESurface(epsilon_n=0.75)
        self.assertAlmostEqual(s10.eps, 0.71625)

        s11 = t.LowESurface(epsilon_n=0.845)
        self.assertAlmostEqual(s11.eps, 0.798525)

    def test_glass_layer(self):

        gl1 = t.GlassLayer(gus=[t.GlassUnit(d=0.003)])
        self.assertAlmostEqual(gl1.r, 0.003)
        self.assertAlmostEqual(gl1.sff.eps, 0.837)
        self.assertAlmostEqual(gl1.sfb.eps, 0.837)

        gl2 = t.GlassLayer(
            gus=[t.GlassUnit(d=0.003)],
            sff=t.LowESurface(epsilon_n=0.15),
            sfb=t.LowESurface(epsilon_n=0.15)
        )
        self.assertAlmostEqual(gl2.sff.eps, 0.168)
        self.assertAlmostEqual(gl2.sfb.eps, 0.168)

        gl3 = t.GlassLayer(gus=[t.GlassUnit(d=0.003), t.GlassUnit(d=0.003)])
        self.assertAlmostEqual(gl3.r, 0.006)

    def test_glass(self):

        g1 = t.Glass(
            gls=[t.GlassLayer(gus=[t.GlassUnit(d=0.003)])],
            als=[]
        )
        self.assertAlmostEqual(g1.r[0], 1 / 20.4013)
        self.assertAlmostEqual(g1.r[1], 0.003)
        self.assertAlmostEqual(g1.r[2], 1 / 8.6198)
        self.assertAlmostEqual(g1.get_h_e(), 20.4013)
        self.assertAlmostEqual(g1.get_h_i(), 8.6198)
        self.assertAlmostEqual(g1.u, 1 / (1 / 20.4013 + 0.003 + 1 / 8.6198))

        ts11, rs11, q11 = g1.get_temp_and_r()
        self.assertAlmostEqual(ts11[0], 5.8343075)
        self.assertAlmostEqual(ts11[1], 6.19138987)
        self.assertAlmostEqual(rs11[0], 0.04901648)
        self.assertAlmostEqual(rs11[1], 0.003)
        self.assertAlmostEqual(rs11[2], 0.11601197)
        self.assertAlmostEqual(q11, 0.0)

        ts12, rs12, q12 = g1.get_temp_and_r(decision_air_layer_temp='calc')
        self.assertAlmostEqual(ts12[0], 5.8343074991436135)
        self.assertAlmostEqual(ts12[1], 6.19138987189045)
        self.assertAlmostEqual(rs12[0], 0.04901648424365114)
        self.assertAlmostEqual(rs12[1], 0.003)
        self.assertAlmostEqual(rs12[2], 0.11601197243555535)
        self.assertAlmostEqual(q12, 0.0)

        ts13, rs13, q13 = g1.get_temp_and_r(decision_air_layer_temp='calc', surface_method='JIS_A2103')
        self.assertAlmostEqual(ts13[0], 4.930922910138927)
        self.assertAlmostEqual(ts13[1], 5.285576406974733)
        self.assertAlmostEqual(rs13[0], 0.0417104832248979)
        self.assertAlmostEqual(rs13[1], 0.003)
        self.assertAlmostEqual(rs13[2], 0.12446873123466247)
        self.assertAlmostEqual(q13, 0.0)

        g2 = t.Glass(
            gls=[
                t.GlassLayer(gus=[t.GlassUnit(d=0.003)]),
                t.GlassLayer(gus=[t.GlassUnit(d=0.003)])
            ],
            als=[
                t.AirLayer(
                    air_property=t.MixedAirProperty(c_air=100.0, c_argon=0.0, c_sf6=0.0, c_krypton=0.0),
                    direction=t.GlassDirection.VERTICAL,
                    s=0.007
                )
            ]
        )

        ts21, rs21, q21 = g2.get_temp_and_r()
        self.assertAlmostEqual(ts21[0], 3.175647357828053)
        self.assertAlmostEqual(ts21[1], 3.370009361151825)
        self.assertAlmostEqual(ts21[2], 12.289531539304013)
        self.assertAlmostEqual(ts21[3], 12.483893542627785)
        self.assertAlmostEqual(rs21[0], 0.04901648424365114)
        self.assertAlmostEqual(rs21[1], 0.003)
        self.assertAlmostEqual(rs21[2], 0.13767385639610663)
        self.assertAlmostEqual(rs21[3], 0.003)
        self.assertAlmostEqual(rs21[4], 0.11601197)
        self.assertAlmostEqual(q21, 0.0)

        ts22, rs22, q22 = g2.get_temp_and_r(decision_air_layer_temp='calc')
        self.assertAlmostEqual(ts22[0], 3.155786765808591)
        self.assertAlmostEqual(ts22[1], 3.3489332234444635)
        self.assertAlmostEqual(ts22[2], 12.337753035938128)
        self.assertAlmostEqual(ts22[3], 12.530899493574)
        self.assertAlmostEqual(rs22[0], 0.04901648424365114)
        self.assertAlmostEqual(rs22[1], 0.003)
        self.assertAlmostEqual(rs22[2], 0.1396166399713072)
        self.assertAlmostEqual(rs22[3], 0.003)
        self.assertAlmostEqual(rs22[4], 0.11601197)
        self.assertAlmostEqual(q22, 0.0)

        ts23, rs23, q23 = g2.get_temp_and_r(decision_air_layer_temp='calc', surface_method='JIS_A2103')
        self.assertAlmostEqual(ts23[0], 2.6973848458545606)
        self.assertAlmostEqual(ts23[1], 2.891001359287642)
        self.assertAlmostEqual(ts23[2], 11.929054072139444)
        self.assertAlmostEqual(ts23[3], 12.122670585572525)
        self.assertAlmostEqual(rs23[0], 0.0417947539395028)
        self.assertAlmostEqual(rs23[1], 0.003)
        self.assertAlmostEqual(rs23[2], 0.14004052473525078)
        self.assertAlmostEqual(rs23[3], 0.003)
        self.assertAlmostEqual(rs23[4], 0.1220556440370072)
        self.assertAlmostEqual(q23, 0.0)

        g3 = t.Glass(
            gls=[
                t.GlassLayer(gus=[t.GlassUnit(d=0.003)]),
                t.GlassLayer(gus=[t.GlassUnit(d=0.003)]),
                t.GlassLayer(gus=[t.GlassUnit(d=0.003)])
            ],
            als=[
                t.AirLayer(
                    air_property=t.MixedAirProperty(c_air=100.0, c_argon=0.0, c_sf6=0.0, c_krypton=0.0),
                    direction=t.GlassDirection.VERTICAL,
                    s=0.007
                ),
                t.AirLayer(
                    air_property=t.MixedAirProperty(c_air=100.0, c_argon=0.0, c_sf6=0.0, c_krypton=0.0),
                    direction=t.GlassDirection.VERTICAL,
                    s=0.007
                )
            ]
        )
        self.assertAlmostEqual(g3.u, 2.2253071434033838)

    def test_glass_with_solar_abs(self):

        g = t.Glass(
            gls=[
                t.GlassLayer(
                    gus=[
                        t.GlassUnit(d=0.003, lmd=1.0),
                        t.GlassUnit(d=0.006, lmd=0.5)
                    ],
                    sff=t.NonLowESurface(),
                    sfb=t.NonLowESurface()
                ),
                t.GlassLayer(
                    gus=[
                        t.GlassUnit(d=0.003, lmd=1.0)
                    ],
                    sff=t.NonLowESurface(),
                    sfb=t.NonLowESurface()
                ),
                t.GlassLayer(
                    gus=[
                        t.GlassUnit(d=0.003, lmd=1.0)
                    ],
                    sff=t.NonLowESurface(),
                    sfb=t.NonLowESurface()
                )
            ],
            als=[
                t.AirLayer(
                    air_property=t.MixedAirProperty(c_air=100.0, c_argon=0.0, c_sf6=0.0, c_krypton=0.0),
                    direction=t.GlassDirection.VERTICAL,
                    s=0.012
                ),
                t.AirLayer(
                    air_property=t.MixedAirProperty(c_air=100.0, c_argon=0.0, c_sf6=0.0, c_krypton=0.0),
                    direction=t.GlassDirection.VERTICAL,
                    s=0.012
                )
            ]
        )

        theta, reg, r_qin = g.get_temp_and_r(
            theta_e=30.0,
            theta_i=25.0,
            surface_method='JIS_A2103',
            season='summer',
            decision_air_layer_temp='calc',
            ia=np.array([9.55935027, 6.8267886, 4.76774099])
        )

        self.assertAlmostEqual(theta[0], 30.26096248)
        self.assertAlmostEqual(theta[1], 30.24131152)
        self.assertAlmostEqual(theta[2], 29.33911229)
        self.assertAlmostEqual(theta[3], 29.31060289)
        self.assertAlmostEqual(theta[4], 27.37569226)
        self.assertAlmostEqual(theta[5], 27.32979106)

        self.assertAlmostEqual(reg[0], 0.07521376)
        self.assertAlmostEqual(reg[1], 0.015)
        self.assertAlmostEqual(reg[2], 0.14815071)
        self.assertAlmostEqual(reg[3], 0.003)
        self.assertAlmostEqual(reg[4], 0.14980114)
        self.assertAlmostEqual(reg[5], 0.003)
        self.assertAlmostEqual(reg[6], 0.1317437)

        self.assertAlmostEqual(r_qin, 8.176926528527648)

if __name__ == '__main__':
    unittest.main()
