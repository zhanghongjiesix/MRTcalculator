from meanradianttemperature import Meanradiant
import pandas as pd
import numpy as np
'''以下路径需要改'''
df = pd.read_excel('E:/博士文件/学习相关/博士论文相关/大论文撰写/图/第六章改/纵横比对太阳辐射传输的影响/辐射数据.xlsx')
multiple_left = []
multiple_right = []
multiple_road = []

for i in np.arange(0, df.shape[0]):
    print('开始计算第{}条的数据'.format(i))
    azimuth = float(df.loc[i, ['azimuth']])
    zenith = float(df.loc[i, ['zenith']])
    solar_beam_energy = float(df.loc[i, ['direct']])
    solar_dif_energy = float(df.loc[i, ['diffuse']])
    a = Radiation_tree(
        canyon_w=40.,
        canyon_h=30.,
        canyon_azimuth=135,
        left_layer=[5., 5., 5., 5., 5., 5.],
        right_layer=[5., 5., 5., 5., 5., 5.],
        w_layer=[5., 5., 5., 5., 5., 5., 5., 5.],
        all_albedo_layer=[0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
                          0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                          0.3, 0.3, 0.3, 0.3, 0.3, 0.3],
        azimuth_s=azimuth, zenith_s=zenith, dir_solar=solar_beam_energy, dif_solar=solar_dif_energy,
        tree=(((10, 4.85), 2.35, 3.9, 1), ((20, 4.85), 2.35, 3.9, 1), ((30, 4.85), 2.35, 3.9, 1)),
    )

    solar_potential_multiple = a.solar_potential()
    left_mean = np.mean(solar_potential_multiple[0: 6])
    road_mean = np.mean(solar_potential_multiple[6: 14])
    right_mean = np.mean(solar_potential_multiple[14: 20])
    multiple_left.append(left_mean)
    multiple_road.append(road_mean)
    multiple_right.append(right_mean)

df.insert(loc=len(df.columns), column='mean_left', value=multiple_left)
df.insert(loc=len(df.columns), column='mean_road', value=multiple_road)
df.insert(loc=len(df.columns), column='mean_right', value=multiple_right)

'''以下路径需要改'''
df.to_excel('E:/博士文件/学习相关/博士论文相关/大论文撰写/图/第六章改/树木种植方式对辐射传输的影响/3行1.xlsx')



