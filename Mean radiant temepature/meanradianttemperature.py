import sys
import numpy as np
import projected_area_factor as paf
import shade_people as sp
import view_factor as vf
sys.path.append('E:\\Python document\\Urban_canyon_model\\Urban_canyon_distribution_with_tree\\Modules')
import tree_radiation

class Meanradiant:
    def __init__(self, canyon_w: float, canyon_w_layer: int, canyon_h: float, canyon_h_layer: int, canyon_azimuth: float, all_albedo_layer: list,
                 azimuth_s: float, zenith_s: float, dir_solar: float, dif_solar: float, all_temp: list,
                 people: list, tree=(((0.3, 10), 0, 0, 0),)):
        self.w = float(canyon_w)  # 街道峡谷宽度
        self.w_layer = int(canyon_w_layer) # 街道峡谷宽度分层
        self.h_layer = int(canyon_h_layer)  # 街道峡谷高度分层
        self.h = float(canyon_h)  # 街道峡谷高度
        self.albedo_w1, self.albedo_g, self.albedo_w2 = all_albedo_layer[0], all_albedo_layer[1], all_albedo_layer[
            2]  # 街道峡谷左墙面、路面、右墙面的反射率
        self.temp_a, self.temp_w1, self.temp_g, self.temp_w2 = all_temp[0], all_temp[1], all_temp[2], all_temp[3]  # 大气、街道峡谷左墙面、路面、右墙面的表面温度
        self.azimuth_c = (canyon_azimuth * np.pi) / 180  # 街道峡谷方位角
        self.azimuth_s = (azimuth_s * np.pi) / 180  # 太阳方位角
        self.zenith_s = (zenith_s * np.pi) / 180  # 太阳天顶角
        self.s_dir = dir_solar  # 直射太阳辐射量
        self.s_dif = dif_solar  # 散射太阳辐射量
        self.x_p, self.y_p = people[0], people[1]  # 人体中心点位置和人体姿态
        self.tree = tree
        assert 0 <= self.azimuth_c < np.pi, '街道峡谷方位角应该在0-pi之间'
        assert 0 <= self.zenith_s < 0.5 * np.pi, '太阳天顶角应该在0-0.5pi之间'
        assert 0 <= self.azimuth_s <= 2 * np.pi, '太阳方位角应在0-2pi之间'
        assert self.azimuth_c != self.azimuth_s, '请避免街道方位角与太阳方位角相同'
        assert 0 < self.albedo_w1 < 1, '表面反射率应该在0-1之间'
        assert 0 < self.albedo_g < 1, '表面反射率应该在0-1之间'
        assert 0 < self.albedo_w2 < 1, '表面反射率应该在0-1之间'
        for i in self.tree:
            assert self.y_p < i[0][1] - i[1], '树冠底部高度应大于人体高度'
            assert self.x_p != i[0][0] + i[1], '请避免人体与树冠的最右侧相切'
            assert self.x_p != i[0][0] - i[1], '请避免人体与树冠的最左侧相切'

    def view_factor(self):
        # 先求没有树木遮挡时人体与街道峡谷各表面之间的视角因子
        view_p_sky, view_p_wall1, view_p_ground, view_p_wall2 = vf.view_wall(w=self.w, h=self.h, people=[self.x_p, self.y_p])
        # 求最近的树木与人体之间的视角因子
        # step1:先找离人体最近的树木：
        distance = []
        for i in self.tree:
            distance_i = ((self.x_p - i[0][0]) ** 2 + (self.y_p - i[0][1]) ** 2) ** 0.5
            distance.append(distance_i)
        a = distance.index(min(distance))
        tree = self.tree[a]
        # step2:求最近树木与人体之间的角系数
        view_p_tree, view_p_treewall1, view_p_treesky, view_p_treewall2 = \
            vf.view_tree(w=self.w, h=self.h, tree=tree, people=[self.x_p, self.y_p])
        # step3:求最近的树的透射率，用来后续求解辐射量
        trans_t = 32 / (9 * tree[2] ** 2) * (1 - np.exp(-0.75 * tree[2])) - 8 / (3 * tree[2]) * np.exp(-0.75 * tree[2])

        return view_p_sky, view_p_wall1, view_p_ground, view_p_wall2, view_p_tree, \
               view_p_treewall1, view_p_treesky, view_p_treewall2, trans_t

    '''
    没有树木的情况下的平均辐射温度的计算
    '''
    def dir_p_notree(self):
        sunfactor = sp.sunfactor(w=self.w, h=self.h, azimuth_canyon=self.azimuth_c, azimuth_solar=self.azimuth_s,
                                 zenith_solar=self.zenith_s, people=[self.x_p, self.y_p])
        # project_factor = paf.paf(self.azimuth_s, self.zenith_s, self.position)
        zenith = (self.zenith_s * 180) / np.pi
        # print(zenith)
        latitude = 90 - zenith
        project_factor = 3.01 * 10 ** (-7) * latitude ** 3 - 6.46 * 10 ** (-5) * latitude ** 2 + 8.34 * 10 ** (-4) * latitude + 0.298
        # print(project_factor)
        dir_people = self.s_dir * project_factor * sunfactor / np.cos(self.zenith_s)
        print('没树遮挡时人体接收到的太阳直射辐射为{}'.format(dir_people))
        return dir_people

    def dif_p_notree(self):
        """
        Returns 没有树木的街道峡谷中人体接收到的太阳散射辐射量
        -------
        """
        # 太阳散射辐射量包括：街道峡谷顶部获得、墙面反射的、路面反射的太阳辐射量
        # step1：使用街道峡谷辐射传输模型计算左墙面、路面和右墙面反射的太阳辐射量
        azimuth_c = (self.azimuth_c * 180) / np.pi
        azimuth_s = (self.azimuth_s * 180) / np.pi
        zenith_s = (self.zenith_s * 180) / np.pi
        r = tree_radiation.Radiation_tree(
            canyon_w=self.w, canyon_h=self.h, canyon_azimuth=azimuth_c,
            left_layer=[self.h / self.h_layer] * self.h_layer,
            w_layer=[self.w / self.w_layer] * self.w_layer,
            right_layer=[self.h / self.h_layer] * self.h_layer,
            all_albedo_layer=[self.albedo_w1] * self.h_layer +
                             [self.albedo_g] * self.w_layer +
                             [self.albedo_w2] * self.h_layer,
            azimuth_s=azimuth_s, zenith_s=zenith_s, dir_solar=self.s_dir, dif_solar=self.s_dif,
        )
        # 到达街道峡谷各个面的太阳辐射量
        solar_wall1, solar_g, solar_wall2 = r.solar_potential()[0], r.solar_potential()[1], r.solar_potential()[2]
        # print('没树时到达的散射辐射量为{}；{}；{}'.format(solar_wall1, solar_g, solar_wall2))
        reflect_solar_wall1 = solar_wall1 * self.albedo_w1
        reflect_solar_ground = solar_g * self.albedo_g
        reflect_solar_wall2 = solar_wall2 * self.albedo_w2
        # 求人体接收到的短波辐射量，首先求视角因子
        vf_sky, vf_wall1, vf_ground, vf_wall2 = vf.view_wall(w=self.w, h=self.h, people=[self.x_p, self.y_p])
        # print('人体与街道峡谷表面的天空视角因子为{}；{}；{}；{}'.format(vf_sky, vf_wall1, vf_ground, vf_wall2))
        dif_sky_p = self.s_dif * vf_sky
        dif_wall1_p = reflect_solar_wall1 * vf_wall1
        dif_ground_p = reflect_solar_ground * vf_ground
        dif_wall2_p = reflect_solar_wall2 * vf_wall2
        print('没树遮挡时人体接收到的散射辐射量分别为{}'.format(dif_sky_p + dif_wall1_p + dif_ground_p + dif_wall2_p))
        return dif_sky_p, dif_wall1_p, dif_ground_p, dif_wall2_p

    def long_p_notree(self):
        """
        Returns 人体接收到的太阳散射辐射量
        -------
        """
        vf_sky, vf_wall1, vf_ground, vf_wall2 = vf.view_wall(w=self.w, h=self.h, people=[self.x_p, self.y_p])
        long_sky_p = 0.95 * 5.67 * 10 ** (-8) * self.temp_a ** 4 * vf_sky
        long_wall1_p = 0.95 * 5.67 * 10 ** (-8) * self.temp_w1 ** 4 * vf_wall1
        long_ground_p = 0.95 * 5.67 * 10 ** (-8) * self.temp_g ** 4 * vf_ground
        long_wall2_p = 0.95 * 5.67 * 10 ** (-8) * self.temp_w2 ** 4 * vf_wall2
        print('没树遮挡时人体接收到的长波辐射量为{}'.format(long_sky_p + long_wall1_p + long_ground_p + long_wall2_p))
        return long_sky_p, long_wall1_p, long_ground_p, long_wall2_p

    def mrt_P_notree(self):
        dir_p = self.dir_p_notree()
        dif_sky_p, dif_wall1_p, dif_ground_p, dif_wall2_p = \
            self.dif_p_notree()[0], self.dif_p_notree()[1], self.dif_p_notree()[2], self.dif_p_notree()[3]
        long_sky_p, long_wall1_p, long_ground_p, long_wall2_p = \
            self.long_p_notree()[0], self.long_p_notree()[1], self.long_p_notree()[2], self.long_p_notree()[3]
        all_radiation_p = (dir_p + dif_sky_p + dif_wall1_p + dif_ground_p + dif_wall2_p) * 0.7 + (long_sky_p + long_wall1_p + long_ground_p + long_wall2_p) * 0.97
        mrt = (all_radiation_p * 10 ** 8 / 0.97 / 5.67) ** 0.25 - 273.15
        print('没树木遮挡时人体mrt为{}'.format(mrt))
        return mrt

    '''
    有树木的情况下的平均辐射温度的计算
    '''

    # 先求人体接收到的太阳直射辐射
    def dir_p_tree(self, treefactor=1):
        """
        Returns 有树木的街道峡谷人体接收到的太阳直射辐射量
        -------
        """
        sunfactor = sp.sunfactor(w=self.w, h=self.h, azimuth_canyon=self.azimuth_c, azimuth_solar=self.azimuth_s,
                                 zenith_solar=self.zenith_s, people=[self.x_p, self.y_p])
        for i in self.tree:
            treefactor_i = sp.treefactor_3d(azimuth_solar=self.azimuth_s, zenith_solar=self.zenith_s, azimuth_canyon=self.azimuth_c,
                                            tree=i, people=[self.x_p, self.y_p])
            treefactor *= treefactor_i
        # project_factor = paf.paf(self.azimuth_s, self.zenith_s, self.position)
        zenith = (self.zenith_s * 180) / np.pi
        # print(zenith)
        latitude = 90 - zenith
        project_factor = 3.01 * 10 ** (-7) * latitude ** 3 - 6.46 * 10 ** (-5) * latitude ** 2 + 8.34 * 10 ** (
            -4) * latitude + 0.298
        # print(project_factor)
        # 算人体得到的辐射量
        dir_people = self.s_dir * project_factor * sunfactor * treefactor / np.cos(self.zenith_s)
        print('有树遮挡时人体接收到的直射辐射：{}'.format(dir_people))
        return dir_people

    def dif_p_tree(self):
        """
        Returns 有树木的街道峡谷中人体接收到的太阳散射辐射量
        -------
        """
        # 太阳散射辐射量包括：街道峡谷顶部获得、墙面反射的、路面反射的太阳辐射量，其中由于树木的遮挡导致人体接收到的太阳散射降低
        # step1：使用街道峡谷辐射传输模型计算左墙面、路面和右墙面反射的太阳辐射量
        azimuth_c = (self.azimuth_c * 180) / np.pi
        azimuth_s = (self.azimuth_s * 180) / np.pi
        zenith_s = (self.zenith_s * 180) / np.pi
        r = tree_radiation.Radiation_tree(
            canyon_w=self.w, canyon_h=self.h, canyon_azimuth=azimuth_c,
            left_layer=[self.h / self.h_layer] * self.h_layer,
            w_layer=[self.w / self.w_layer] * self.w_layer,
            right_layer=[self.h / self.h_layer] * self.h_layer,
            all_albedo_layer=[self.albedo_w1] * self.h_layer +
                             [self.albedo_g] * self.w_layer +
                             [self.albedo_w2] * self.h_layer,
            azimuth_s=azimuth_s, zenith_s=zenith_s, dir_solar=self.s_dir, dif_solar=self.s_dif,
            tree=self.tree
        )
        # 求街道峡谷各个面反射的太阳辐射量
        solar_wall1, solar_g, solar_wall2 = r.solarpotential_mean()[0], r.solarpotential_mean()[1], r.solarpotential_mean()[2]
        # print('有树时到达的散射辐射量为{}；{}；{}'.format(solar_wall1, solar_g, solar_wall2))
        reflect_solar_wall1 = solar_wall1 * self.albedo_w1
        reflect_solar_ground = solar_g * self.albedo_g
        reflect_solar_wall2 = solar_wall2 * self.albedo_w2
        view_p_sky, view_p_wall1, view_p_ground, view_p_wall2, view_p_tree, \
        view_p_treewall1, view_p_treesky, view_p_treewall2, trans_t = self.view_factor()
        # 求人体接收到的短波辐射量
        dif_sky_p = self.s_dif * view_p_treesky * trans_t + self.s_dif * (view_p_sky - view_p_treesky) * trans_t
        dif_wall1_p = reflect_solar_wall1 * view_p_treewall1 * trans_t + reflect_solar_wall1 * (view_p_wall1 - view_p_treewall1) * trans_t
        dif_ground_p = reflect_solar_ground * view_p_ground
        dif_wall2_p = reflect_solar_wall2 * view_p_treewall2 * trans_t + reflect_solar_wall2 * (view_p_wall2 - view_p_treewall2) * trans_t
        dif_P = dif_sky_p + dif_wall1_p + dif_wall2_p + dif_ground_p
        print('有树遮挡时人体接收到的散射辐射：{}'.format(dif_P))
        return dif_P

    def long_p_tree(self):
        """
        Returns 有树木的街道峡谷中人体接收到的长波辐射量
        -------
        """
        # step1: 先把视角因子都列出来
        view_p_sky, view_p_wall1, view_p_ground, view_p_wall2, view_p_tree, \
        view_p_treewall1, view_p_treesky, view_p_treewall2, trans_t = self.view_factor()
        long_sky_p = (view_p_sky - view_p_treesky) * 0.95 * 5.67 * 10 ** (-8) * self.temp_a ** 4
        long_tree_p = view_p_tree * 0.95 * 5.67 * 10 ** (-8) * self.temp_a ** 4
        long_wall1_p = (view_p_wall1 - view_p_treewall1) * 0.95 * 5.67 * 10 ** (-8) * self.temp_w1 ** 4
        long_ground_p = view_p_ground * 0.95 * 5.67 * 10 ** (-8) * self.temp_g ** 4
        long_wall2_p = (view_p_wall2 - view_p_treewall2) * 0.95 * 5.67 * 10 ** (-8) * self.temp_w2 ** 4
        long_p = long_sky_p + long_tree_p + long_wall1_p + long_wall2_p + long_ground_p
        print('有树木遮挡时人体接收到的长波辐射量为{}'.format(long_p))
        return long_p

    def mrt_P_tree(self):
        dir_p = self.dir_p_tree()
        dif_p = self.dif_p_tree()
        long_p = self.long_p_tree()
        all_radiation_p = (dir_p + dif_p) * 0.7 + long_p * 0.97
        mrt = (all_radiation_p * 10 ** 8 / 0.97 / 5.67) ** 0.25 - 273.15
        print('有树木遮挡时人体mrt为{}'.format(mrt))
        return mrt

if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    '''以下路径需要改'''
    df = pd.read_excel('E:/博士文件/学习相关/博士论文相关/大论文撰写/图/第六章改/观音桥步行街人体MRT分析/夏季输入数据.xlsx')
    MRT = []

    for i in np.arange(0, df.shape[0]):
        print('开始计算第{}条的数据'.format(i))
        azimuth = float(df.loc[i, ['azimuth']])
        zenith = float(df.loc[i, ['zenith']])
        solar_beam_energy = float(df.loc[i, ['direct']])
        solar_dif_energy = float(df.loc[i, ['diffuse']])
        air, left, wall, right = float(df.loc[i, ['air_temp']]), float(df.loc[i, ['left_temp']]), float(df.loc[i, ['wall_temp']]), float(df.loc[i, ['right_temp']])
        mrt_people = Meanradiant(
                canyon_w=40, canyon_w_layer=8,
                canyon_h=30, canyon_h_layer=6,
                canyon_azimuth=135, all_albedo_layer=[0.3, 0.2, 0.3],
                all_temp=[air, left, wall, right],
                azimuth_s=azimuth, zenith_s=zenith, dif_solar=solar_dif_energy, dir_solar=solar_beam_energy, people=[21, 1.1]
                , tree=(((20, 8.0), 6.0, 4.77, 1), )
            )
        result = mrt_people.mrt_P_notree()
        MRT.append(result)
    df.insert(loc=len(df.columns), column='MRT', value=MRT)
    '''以下路径需要改'''
    df.to_excel('E:/博士文件/学习相关/博士论文相关/大论文撰写/图/第六章改/树木种类对MRT的影响/榕树1.xlsx')

    # mrt_people = Meanradiant(
    #     canyon_w=20, canyon_w_layer=8,
    #     canyon_h=20, canyon_h_layer=8,
    #     canyon_azimuth=0, all_albedo_layer=[0.3, 0.2, 0.3],
    #     all_temp=[311.4, 331.2, 319.9, 318.9],
    #     azimuth_s=96.3, zenith_s=42.9, dif_solar=102.5, dir_solar=476.5, people=[10, 1.1]
    #     # , tree=(((1, 1.15), 0.735, 2.49, 1),)
    # )
    # a = mrt_people.mrt_P_notree()
    # # b = mrt_people.mrt_P_tree()


