import numpy as np
import sys
sys.path.append('E:\\Python document\\Urban_canyon_model\\Urban_canyon_distribution_with_tree\\Modules')
import tree_radiation




# '''''''''''''''''''''''''''''''''''''''''''''先把视角因子等需要的参数求出来'''''''''''''''''''''''''''''''''''''''''''''
r = tree_radiation.Radiation_tree(canyon_w=5, canyon_h=20, canyon_azimuth=0,
                                  left_layer=[2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5],
                                  w_layer=[2.5, 2.5],
                                  right_layer=[2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5],
                                  all_albedo_layer=[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
                                                    0.2, 0.2,
                                                    0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3],
                                  azimuth_s=284.9, zenith_s=81.5, dif_solar=6.9, dir_solar=32.1,
                                  # tree=(((1, 0.75), 0.25, 2.49, 1),)
                                  )
wind = 0.78      # 风速
Ta = 312.3        # 空气温度
error = 0.1    # 误差值
T1 = 300        # 初始wall1温度
Tg = 300        # 初始ground温度
T2 = 300        # 初始wall2温度
# '''''''''''''''''''''''''''''''''''''''''''''''''需要修改的参数'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
views = r.view_mat_temp()
print(r.solar_potential())
solars = r.solarpotential_mean()
view1_g, view1_2, view1_s = views[1], views[2], views[3]
viewg_1, viewg_2, viewg_s = views[4], views[6], views[7]
view2_1, view2_g, view2_s = views[8], views[9], views[11]
views_1, views_g, views_2 = views[12], views[13], views[14]
solar_1 = solars[0] * (1 - 0.3)
solar_g = solars[1] * (1 - 0.2)             # 改反射率
solar_2 = solars[2] * (1 - 0.3)
print('the solar absorbed by wall1, ground and wall2 is {},{},{}'.format(solar_1, solar_g, solar_2))

# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# '''''''''''''''''''''''''''''''''''''''''''''需要的函数'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def F1(T1, Tg, T2, solar_1, wind, a1):
    return Ta + a1 * (0.95 * 5.67 * 10 ** (-8) * (
                Tg ** 4 * viewg_1 + T2 ** 4 * view2_1 + Ta ** 4 * views_1 - T1 ** 4) + solar_1) / (
                       1.05 * 6.2 + 4.26 * wind) - T1
def dF11(T1, wind, a1):
    return (a1 * 0.95 * 5.67 * 10 ** (-8) * (- 4 * T1 ** 3)) / (1.05 * 6.2 + 4.26 * wind) - 1
def dF1g(Tg, viewg_1, wind, a1):
    return (a1 * 0.95 * 5.67 * 10 ** (-8) * 4 * Tg ** 3 * viewg_1) / (1.05 * 6.2 + 4.26 * wind)
def dF12(T2, view2_1, wind, a1):
    return (a1 * 0.95 * 5.67 * 10 ** (-8) * 4 * T2 ** 3 * view2_1) / (1.05 * 6.2 + 4.26 * wind)
def F2(T1, Tg, T2, solar_2, wind, a2):
    return Ta + a2 * (0.95 * 5.67 * 10 ** (-8) * (
                T1 ** 4 * view1_2 + Tg ** 4 * viewg_2 + Ta ** 4 * views_2 - T2 ** 4) + solar_2) / (
                       1.05 * 6.2 + 4.26 * wind) - T2
def dF22(T2, wind, a2):
    return (a2 * 0.95 * 5.67 * 10 ** (-8) * (- 4 * T2 ** 3)) / (1.05 * 6.2 + 4.26 * wind) - 1
def dF2g(Tg, viewg_2, wind, a2):
    return (a2 * 0.95 * 5.67 * 10 ** (-8) * 4 * Tg ** 3 * viewg_2) / (1.05 * 6.2 + 4.26 * wind)
def dF21(T1, view1_2, wind, a2):
    return (a2 * 0.95 * 5.67 * 10 ** (-8) * 4 * T1 ** 3 * view1_2) / (1.05 * 6.2 + 4.26 * wind)
def Fg(T1, Tg, T2, solar_g, wind, ag):
    return Ta + ag * (0.95 * 5.67 * 10 ** (-8) * (
                T1 ** 4 * view1_g + T2 ** 4 * view2_g + Ta ** 4 * views_g - Tg ** 4) + solar_g) / (
                       1.05 * 6.2 + 4.26 * wind) - Tg


def dFgg(Tg, wind, ag):
    return (ag * 0.95 * 5.67 * 10 ** (-8) * (- 4 * Tg ** 3)) / (1.05 * 6.2 + 4.26 * wind) - 1
def dFg1(T1, view1_g, wind, ag):
    return (ag * 0.95 * 5.67 * 10 ** (-8) * 4 * T1 ** 3 * view1_g) / (1.05 * 6.2 + 4.26 * wind)
def dFg2(T2, view2_g, wind, ag):
    return (ag * 0.95 * 5.67 * 10 ** (-8) * 4 * T2 ** 3 * view2_g) / (1.05 * 6.2 + 4.26 * wind)


# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# ''''''''''''''''''''''''''''''''''''''''''''''判断系数因子用的'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
factor1 = 0.95 * 5.67 * 10 ** (-8) * (Tg ** 4 * viewg_1 + T2 ** 4 * view2_1 + Ta ** 4 * views_1 - T1 ** 4) + solar_1
factor2 = 0.95 * 5.67 * 10 ** (-8) * (T1 ** 4 * view1_2 + Tg ** 4 * viewg_2 + Ta ** 4 * views_2 - T2 ** 4) + solar_2
factorg = 0.95 * 5.67 * 10 ** (-8) * (T1 ** 4 * view1_g + T2 ** 4 * view2_g + Ta ** 4 * views_g - Tg ** 4) + solar_g
if factor1 > 0:
    a1 = 0.81
else:
    a1 = 0.68
if factor2 > 0:
    a2 = 0.81
else:
    a2 = 0.68
if factorg > 0:
    ag = 0.81
else:
    ag = 0.68

ini_T = np.array([[T1], [Tg], [T2]])
b = np.array([F1(T1, Tg, T2, solar_1, wind, a1) ** 2, F2(T1, Tg, T2, solar_2, wind, a2) ** 2,
              Fg(T1, Tg, T2, solar_g, wind, ag) ** 2])
while b.all() > error:
    arr_1 = np.array(
        [[F1(T1, Tg, T2, solar_1, wind, a1)], [F2(T1, Tg, T2, solar_2, wind, a2)], [Fg(T1, Tg, T2, solar_g, wind, ag)]])
    arr_2 = np.linalg.inv(np.array([[dF11(T1, wind, a1), dF12(T2, view2_1, wind, a1), dF1g(Tg, viewg_1, wind, a1)],
                                    [dF21(T1, view1_2, wind, a2), dF22(T2, wind, a2), dF2g(Tg, viewg_2, wind, a2)],
                                    [dFg1(T1, view1_g, wind, ag), dFg2(T2, view2_g, wind, ag), dFgg(Tg, wind, ag)]]))
    ini_T = ini_T - np.dot(arr_2, arr_1)
    # print('newT{}'.format(ini_T))
    T1 = ini_T[0][0]
    T2 = ini_T[1][0]
    Tg = ini_T[2][0]
    b = np.array([F1(T1, Tg, T2, solar_1, wind, a1) ** 2, F2(T1, Tg, T2, solar_2, wind, a2) ** 2,
                  Fg(T1, Tg, T2, solar_g, wind, ag) ** 2])
    factor1 = 0.95 * 5.67 * 10 ** (-8) * (Tg ** 4 * viewg_1 + T2 ** 4 * view2_1 + Ta ** 4 * views_1 - T1 ** 4) + solar_1
    factor2 = 0.95 * 5.67 * 10 ** (-8) * (T1 ** 4 * view1_2 + Tg ** 4 * viewg_2 + Ta ** 4 * views_2 - T2 ** 4) + solar_2
    factorg = 0.95 * 5.67 * 10 ** (-8) * (T1 ** 4 * view1_g + T2 ** 4 * view2_g + Ta ** 4 * views_g - Tg ** 4) + solar_g
    if factor1 > 0:
        a1 = 0.81
    else:
        a1 = 0.68
    if factor2 > 0:
        a2 = 0.81
    else:
        a2 = 0.68
    if factorg > 0:
        ag = 0.81
    else:
        ag = 0.68
print('the final temperature of wall1, ground and wall2 is {},{},{}'.format(T1, Tg, T2))


