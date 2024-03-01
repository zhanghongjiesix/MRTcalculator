import numpy as np


def sunfactor(w, h, azimuth_canyon, azimuth_solar, zenith_solar, people):
    assert 0 <= azimuth_canyon <= np.pi, '街道峡谷的方位角应该在0~pi之间'
    assert 0 <= azimuth_solar <= 2 * np.pi, '太阳方位角应该在0~2pi之间'
    assert 0 <= zenith_solar <= np.pi, '太阳天顶角应该在0~pi之间'
    xp, yp = people[0], people[1]
    # 先求光线的当量天顶角：
    angle = np.arctan(np.tan(zenith_solar) * np.sin(azimuth_solar - azimuth_canyon))
    slope = np.tan(np.pi / 2 - angle)
    if azimuth_canyon < azimuth_solar < azimuth_canyon + np.pi:  # 此时太阳位于街道峡谷的右侧，右墙面顶点决定阴影位置
        # 求此时右墙面顶点形成的阴影的公式:y = k (x - w) + h，
        y = slope * (xp - w) + h  # 求当x等于xp时对应的y值，如果y值大于yp，那么行人位于阴影之中，如果y值小于yp，那么行人位于阳光之中
    else:  # 此时太阳位于街道峡谷的左侧，左墙面顶点决定阴影位置
        # 求此时左墙面顶点形成的阴影的公式：y = kx + h
        y = slope * xp + h  # 求当x等于xp时对应的y值，如果y值大于yp，那么行人位于阴影之中，如果y值小于yp，那么行人位于阳光之中
    if y > yp:
        return 0
    else:
        return 1


def treefactor_2d(azimuth_canyon, azimuth_solar, zenith_solar, tree, people):
    assert 0 <= azimuth_canyon <= np.pi, '街道峡谷的方位角应该在0~pi之间'
    assert 0 <= azimuth_solar <= 2 * np.pi, '太阳方位角应该在0~2pi之间'
    assert 0 <= zenith_solar <= np.pi, '太阳天顶角应该在0~pi之间'
    xp, yp = people[0], people[1]
    xt, yt, r, LAI = tree[0][0], tree[0][1], tree[1], tree[2]
    # 判断人体重心点是不是在树冠的阴影之中，使用人体到树冠阴影上下切线的距离进行判断，如果点到两个直线的距离等树冠的直径，那么人体就位于树冠阴影之中。
    # 先求光线的当量天顶角：
    angle = np.arctan(np.tan(zenith_solar) * np.sin(azimuth_solar - azimuth_canyon))
    slope = np.tan(np.pi / 2 - angle)
    if azimuth_canyon < azimuth_solar < azimuth_canyon + np.pi:
        # 根据树冠的中线点（xt,yt）和太阳光线的斜率k求过树冠中线点的直线方程：y-yt = k(x-xt)
        # 上切线的直线方程为：y = k(x - xt + abs(r * cos(angle))) + yt + abs(r * sin(angle))
        # 下切线的直线方程为：y = k(x - xt - abs(r * cos(angle))) + yt - abs(r * sin(angle))
        # '点到上切线的距离由以下公式计算：'
        d1 = abs(slope * xp - yp + slope * abs(r * np.cos(angle)) + yt + abs(r * np.sin(angle)) - slope * xt) / (
                    slope ** 2 + 1) ** 0.5
        # print(d1)
        # '点到下切线的距离由以下公式计算：'
        d2 = abs(slope * xp - yp - slope * abs(r * np.cos(angle)) + yt - abs(r * np.sin(angle)) - slope * xt) / (
                    slope ** 2 + 1) ** 0.5
        # print(d2)
    else:
        # 上切线的直线方程为：y = k(x - xt - r * abs(cos(angle))) + yt + r * abs(sin(angle))
        # 下切线的直线方程为：y = k(x - xt + r * abs(cos(angle))) + yt - r * abs(sin(angle))
        # 点到上切线的距离由以下公式计算
        d1 = abs(slope * xp - yp - slope * abs(r * np.cos(angle)) + yt + abs(r * np.sin(angle)) - slope * xt) / (
                    slope ** 2 + 1) ** 0.5
        # print(d1)
        # 点到下切线的距离由以下公式计算
        d2 = abs(slope * xp - yp + slope * abs(r * np.cos(angle)) + yt - abs(r * np.sin(angle)) - slope * xt) / (
                    slope ** 2 + 1) ** 0.5
        # print(d2)
    if abs(d1 + d2 - 2 * r) <= 0.05:
        return 32 / (9 * LAI ** 2) * (1 - np.exp(-0.75 * LAI)) - 8 / (3 * LAI) * np.exp(-0.75 * LAI)
    else:
        return 1

def treefactor_3d(azimuth_canyon, azimuth_solar, zenith_solar, tree, people):
    assert 0 <= azimuth_canyon <= np.pi, '街道峡谷的方位角应该在0~pi之间'
    assert 0 <= azimuth_solar <= 2 * np.pi, '太阳方位角应该在0~2pi之间'
    assert 0 <= zenith_solar <= 0.5 * np.pi, '太阳天顶角应该在0~0.5pi之间'
    xp, yp, zp = people[0], 0, people[1]
    xt, yt, zt, r, LAI = tree[0][0], 0, tree[0][1], tree[1], tree[2]
    # 求光线的方向向量
    l = np.sin(zenith_solar) * np.sin(azimuth_solar - azimuth_canyon)
    m = np.sin(zenith_solar) * np.cos(azimuth_solar - azimuth_canyon)
    n = np.cos(zenith_solar)
    # 判断人体是否在树冠阴影之中
    distance_1 = (xp - xt) ** 2 + (yp - yt) ** 2 + (zp - zt) ** 2 - r ** 2
    distance_2 = (l * (xp - xt) + m * (yp - yt) + n * (zp - zt)) ** 2 / (l ** 2 + m ** 2 + n ** 2)
    if distance_1 <= distance_2:
        return 32 / (9 * LAI ** 2) * (1 - np.exp(-0.75 * LAI)) - 8 / (3 * LAI) * np.exp(-0.75 * LAI)
    else:
        return 1



if __name__ == '__main__':
    azimuth_solar = 278.6 * np.pi / 180
    zenith_solar = 86.3 * np.pi / 180
    a = treefactor_2d(azimuth_canyon=0, azimuth_solar=azimuth_solar, zenith_solar=zenith_solar, tree=((4, 6), 3, 3, 0),
                      people=[5, 1])
    print(a)
    b = treefactor_3d(azimuth_canyon=0, azimuth_solar=azimuth_solar, zenith_solar=zenith_solar, tree=((4, 6), 3, 3, 0),
                      people=[5, 1])
    print(b)