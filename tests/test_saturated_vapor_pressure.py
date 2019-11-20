import unittest

from src.saturated_vapor_pressure import get_saturated_vapor_pressure as svp
from src.saturated_vapor_pressure import get_saturated_vapor_pressure_differential as dsvp_dt


class TestSVP(unittest.TestCase):
    
    def test_get_svp(self):
        
        test_patterns = [
                ('SONNTAG',      'water', 273.15,    611.212840000,   44.409707360),
                ('SONNTAG',      'water', 323.15,  12352.743250000,  612.976269100),
                ('SONNTAG',      'water', 373.15, 101419.041700000, 3619.570360000),
                ('SONNTAG',      'ice'  , 223.15,      3.935793929,    0.486049413),
                ('SONNTAG',      'ice'  , 273.15,    611.153849700,   50.323260310),
                ('WMO',          'water', 273.15,    610.695095700,   44.376240950),
                ('WMO',          'water', 323.15,  12338.999600000,  612.272861000),
                ('WMO',          'water', 373.15, 101325.129100000, 3616.814800000),
                ('WMO',          'ice'  , 223.15,      6.354201838,    0.728414468),
                ('WMO',          'ice'  , 273.15,    610.695095700,   44.376240950),
                ('WexlerHyland', 'water', 273.15,    611.212867500,   44.401742580),
                ('WexlerHyland', 'water', 323.15,  12349.856470000,  612.866982700),
                ('WexlerHyland', 'water', 373.15, 101418.716800000, 3619.664421000),
                ('WexlerHyland', 'ice'  , 223.15,      3.938985632,    0.486391125),
                ('WexlerHyland', 'ice'  , 273.15,    611.153570900,   50.326485290),
                ('Tetens',       'water', 273.15,    610.335650900,   44.420613910),
                ('Tetens',       'water', 323.15,  12328.919330000,  612.151168500),
                ('Tetens',       'water', 373.15, 102157.027200000, 3679.899908000),
                ('Tetens',       'ice'  , 223.15,      3.812109171,    0.476776572),
                ('Tetens',       'ice'  , 273.15,    610.276966400,   50.284537900),
                ('BriggsSacket', 'water', 273.15,    613.027145500,   44.476649680),
                ('BriggsSacket', 'water', 323.15,  12373.380390000,  613.621835400),
                ('BriggsSacket', 'water', 373.15, 101447.181100000, 3623.833972000),
                ('BriggsSacket', 'ice'  , 223.15,      3.963611731,    0.490459045),
                ('BriggsSacket', 'ice'  , 273.15,    613.015817400,   50.360785400),       
                ('Antoine',      'water', 273.15,    604.987807200,   44.370885400),
                ('Antoine',      'water', 323.15,  12341.583900000,  612.074421200),
                ('Antoine',      'water', 373.15, 101349.548100000, 3624.104014000),
                ('Antoine',      'ice'  , 223.15,      5.625853212,    0.671408508),
                ('Antoine',      'ice'  , 273.15,    604.987807200,   44.370885400),
                ('GoffGratch',   'water', 273.15,    610.663250200,   44.374468060),
                ('GoffGratch',   'water', 323.15,  12338.741550000,  612.265913200),
                ('GoffGratch',   'water', 373.15, 101325.000000000, 3616.842993000),
                ('GoffGratch',   'ice'  , 223.15,      3.936454864,    0.486146856),
                ('GoffGratch',   'ice'  , 273.15,    611.226430100,   50.338119000),
                ]
        
        for equation, status, t, expected_svp, expected_dsvp_dt in test_patterns:
            with self.subTest(equation=equation, status=status, t=t):
                actual_svp = svp(equation, status, t)
                actual_dsvp_dt = dsvp_dt(equation, status, t)
                self.assertAlmostEqual(actual_svp, expected_svp, delta=expected_svp*0.0001)
                self.assertAlmostEqual(actual_dsvp_dt, expected_dsvp_dt, delta=expected_svp*0.0001)


if __name__ == "__main__":
    unittest.main()
