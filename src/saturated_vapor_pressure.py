import math
import matplotlib.pyplot as plt
from collections import namedtuple

def saturated_vapor_pressure_SONNTAG(status, T):
    Coeff = namedtuple('Coeff',('a1','a2','a3','a4','a5'))
    c = {
        'water' : Coeff( -6096.9385, 21.2409642, -0.02711193,   0.00001673952,   2.433502   ),
        'ice'   : Coeff( -6024.5282, 29.32707,    0.010613863, -0.000013198825, -0.49382577 )
    }[status]
    k = c.a1 / T + c.a2 + c.a3 * T + c.a4 * T**2 + c.a5 * math.log(T)
    pvs = math.exp(k)
    dpvs_dT = pvs * ( - c.a1 / ( T**2 ) + c.a3 + 2 * c.a4 * T + c.a5 / T )
    return (pvs,dpvs_dT)


def saturated_vapor_pressure_WMO(T):
    ew = 2.78614 + 10.79574 * ( 1.0 - 273.16 / T ) \
       - 5.028 * math.log10( T / 273.16 ) \
       + 1.50475 * 10**(-4) * ( 1.0 - 10.0**( -8.2969 * ( T / 273.16 - 1.0 ) ) ) \
       + 0.42873 * 10**(-3) * ( 10**( 4.76955 * ( 1.0 - 273.16 / T ) )  - 1.0 )
    dew_dT = 10.79574 * 273.16 / (T**2) \
           - 5.028 / math.log(10) / T \
           + 1.50475 * 10**(-4) * 8.2969 / 273.16 * math.log(10.0) * 10.0**( -8.2969 * (T/273.16 - 1.0) ) \
           + 0.42873 * 10**(-3) * 4.76955 * 273.16 / (T**2) * math.log(10.0) * 10**( 4.76955 * (1.0-273.16/T) )
    return ( 10**ew, 10**ew * dew_dT * math.log(10.0) )


def saturated_vapor_pressure_WH(status, T):
    Coeff = namedtuple('Coeff', ('a1','b1','a2','b2','a3','b3','a4','b4','a5','b5','a6','b6','a7','b7'))
    c = {
        'water' : Coeff( -0.58002206, 4, 0.13914993, 1, -0.48640239, -1, 0.41764768, -4, -0.14452093, -7,  0.0,          0, 0.65459673, 1),
        'ice'   : Coeff( -0.56745359, 4, 0.63925247, 1, -0.96778430, -2, 0.62215701, -6,  0.20747825, -8, -0.94840240, -12, 0.41635019, 1)
    }[status]
    k = c.a1 * 10**c.b1 / T + c.a2 * 10**c.b2 + c.a3 * 10**c.b3 * T + c.a4 * 10**c.b4 * T**2 + c.a5 * 10**c.b5 * T**3 + c.a6 * 10**c.b6 * T**4 + c.a7 * 10**c.b7 * math.log(T)
    pvs     = math.exp(k)
    dpvs_dT = pvs * ( -c.a1 * 10**(c.b1) / T**2 + c.a3 * 10**(c.b3) + 2 * c.a4 * 10**c.b4 * T + 3 * c.a5 * 10**c.b5 * T**2 + 4 * c.a6 * 10**c.b6 * T**3 + c.a7 * 10**c.b7 / T )
    return ( pvs, dpvs_dT )


def saturated_vapor_pressure_tetens(status, T):
    Coeff = namedtuple('Coeff',('a','b'))
    c = {
        'water' : Coeff( 17.2693882, 35.86 ),
        'ice'   : Coeff( 21.8745584,  7.66 )
    }[status]
    pvs     = 6.1078 * 10**2 * math.exp( c.a * ( T - 273.16 ) / ( T - c.b ) )
    dpvs_dT = pvs * ( c.a / ( T - c.b ) - c.a * ( T - 273.16 ) / (( T - c.b )**2) )
    return ( pvs, dpvs_dT )


def saturated_vapor_pressure_BS(status, T):
    Coeff = namedtuple('Coeff', ('a1','a2','a3','a4','a5'))
    c = {
        'water' : Coeff( -2313.0338, -164.03307,  38.053682, -0.13844344, 0.000074465367),
        'ice'   : Coeff( -5631.1206,   -8.363602,  8.2312,   -0.03861449, 0.0000277494)
    }[status]
    k = c.a1 / T + c.a2 + c.a3 * math.log(T) + c.a4 * T + c.a5 * T**2 - math.log(10.0)
    pvs     = math.exp(k)
    dpvs_dT = pvs * ( -c.a1 / T**2 + c.a3 / T + c.a4 + 2 * c.a5 * T )
    return (pvs, dpvs_dT)



def saturated_vapor_pressure_antoine(T):
    A = 8.02754
    B = 1705.616
    C = 231.405 - 273.15
    pvs = 10**(A - (B/(T+C)))
    dpvs_dt = pvs * math.log(10) * ( B / ( T + C )**2 )
    return (pvs * 101325 / 760, dpvs_dt * 101325 / 760) # 760mmHg = 101325 Pa


def saturated_vapor_pressure_GoffGratch(status, T):
    a_1  =   -7.90298
    a_2  =    5.02808
    a_3  =   -1.3816 * 10**(-7)
    a_4  =   11.344
    a_5  =    8.1328 * 10**(-3)
    a_6  =   -3.49149
    T_st =  373.15
    e_st = 1013.25
    k_w  = a_1 * ( T_st/T - 1 ) + a_2 * math.log10( T_st/T ) + a_3 * ( 10**(a_4*(1-T/T_st)) - 1 ) + a_5 * ( 10**(a_6*(T_st/T-1)) - 1 ) + math.log10(e_st)
    pvs_w = 10**k_w
    dpvs_dT_w = pvs_w * math.log(10) * ( - a_1 * T_st / T**2 - a_2 / T / math.log(10) - a_3 * math.log(10) * 10**(a_4*(1-T/T_st)) * a_4 / T_st - a_5 * math.log(10) * 10**(a_6*(T_st/T-1)) * a_6 * T_st / T**2 )
    b_1  =  -9.09718
    b_2  =  -3.56654
    b_3  =   0.876793
    T_0  = 273.16
    e_i0 =   6.1173
    k_i  = b_1 * ( T_0/T - 1 ) + b_2 * math.log10( T_0/T ) + b_3 * (1-T/T_0) + math.log10(e_i0)
    pvs_i = 10**k_i
    dpvs_dT_i = pvs_i * math.log(10) * ( -b_1 * T_0 / T**2 - b_2 / math.log(10) / T - b_3 / T_0 )
    return {
        'water' : (pvs_w * 100, dpvs_dT_w * 100),
        'ice'   : (pvs_i * 100, dpvs_dT_i * 100)
    }[status]


def get_saturated_vapor_pressure( equation, status, T ):
    if T < 0:
        raise ValueError('ERROR: Temperature can not be less than 0 K.')
    if equation == 'SONNTAG':
        return saturated_vapor_pressure_SONNTAG(status,T)[0]
    elif equation == 'WMO':
        return saturated_vapor_pressure_WMO(T)[0]
    elif equation == 'WexlerHyland':
        return saturated_vapor_pressure_WH(status,T)[0]
    elif equation == 'Tetens':
        return saturated_vapor_pressure_tetens(status,T)[0]
    elif equation == 'BriggsSacket':
        return saturated_vapor_pressure_BS(status,T)[0]
    elif equation == 'Antoine':
        return saturated_vapor_pressure_antoine(T)[0]
    elif equation == 'GoffGratch':
        return saturated_vapor_pressure_GoffGratch(status,T)[0]
    else:
        raise ValueError('ERROR: false name of saturated vapor equation')

#    return {
#        'SONNTAG'      : saturated_vapor_pressure_SONNTAG(status,T)[0],
#        'WMO'          : saturated_vapor_pressure_WMO(T)[0],
#        'WexlerHyland' : saturated_vapor_pressure_WH(status,T)[0],
#        'Tetens'       : saturated_vapor_pressure_tetens(status,T)[0],
#        'BriggsSacket' : saturated_vapor_pressure_BS(status,T)[0],
#        'Antoine'      : saturated_vapor_pressure_antoine(T)[0],
#        'GoffGratch'   : saturated_vapor_pressure_GoffGratch(status,T)[0]
#    }[equation]




def get_saturated_vapor_pressure_differential( equation, status, T ):
    if T < 0:
        raise ValueError('ERROR: Temperature can not be less than 0 K.')
    return {
        'SONNTAG'      : saturated_vapor_pressure_SONNTAG(status,T)[1],
        'WMO'          : saturated_vapor_pressure_WMO(T)[1],
        'WexlerHyland' : saturated_vapor_pressure_WH(status,T)[1],
        'Tetens'       : saturated_vapor_pressure_tetens(status,T)[1],
        'BriggsSacket' : saturated_vapor_pressure_BS(status,T)[1],
        'Antoine'      : saturated_vapor_pressure_antoine(T)[1],
        'GoffGratch'   : saturated_vapor_pressure_GoffGratch(status,T)[1]
    }[equation]