import numpy as np
from sympy import *
import math

def view_wall(w, h, people):
    # 根据w和h得到街道峡谷的四个顶点：从左上到右上分别为：
    x1, y1 = 0, h
    x2, y2 = 0, 0
    x3, y3 = w, 0
    x4, y4 = w, h
    xp, yp = people[0], people[1]
    # 求人体与天空角系数
    angle_p_sky = math.acos(((x1 - xp) ** 2 + (y1 - yp) ** 2 + (x4 - xp) ** 2 + (y4 - yp) ** 2 - w ** 2)/(
                  2 * ((x1 - xp) ** 2 + (y1 - yp) ** 2) ** 0.5 *((x4 - xp) ** 2 + (y4 - yp) ** 2) ** 0.5))
    view_p_sky = angle_p_sky / 2 / np.pi
    # 求人体与左墙面的角系数
    angle_p_wall1 = math.acos(((x1 - xp) ** 2 + (y1 - yp) ** 2 + (x2 - xp) ** 2 + (y2 - yp) ** 2 - h ** 2)/(
                    2 * ((x1 - xp) ** 2 + (y1 - yp) ** 2) ** 0.5 * ((x2 - xp) ** 2 + (y2 - yp) ** 2) ** 0.5))
    view_p_wall1 = angle_p_wall1 / 2 / np.pi
    # 求人体与路面的角系数
    angle_p_ground = math.acos(((x2 - xp) ** 2 + (y2 - yp) ** 2 + (x3 - xp) ** 2 + (y3 - yp) ** 2 - w ** 2)/(
                     2 * ((x2 - xp) ** 2 + (y2 - yp) ** 2) ** 0.5 *((x3 - xp) ** 2 + (y3 - yp) ** 2) ** 0.5))
    view_p_ground = angle_p_ground / 2 / np.pi
    # 求人体与右墙面的角系数
    angle_p_wall2 = math.acos(((x3 - xp) ** 2 + (y3 - yp) ** 2 + (x4 - xp) ** 2 + (y4 - yp) ** 2 - h ** 2) / (
                    2 * ((x3 - xp) ** 2 + (y3 - yp) ** 2) ** 0.5 * ((x4 - xp) ** 2 + (y4 - yp) ** 2) ** 0.5))
    view_p_wall2 = angle_p_wall2 / 2 / np.pi
    return view_p_sky, view_p_wall1, view_p_ground, view_p_wall2

def line(tree, people):
    xp, yp = people[0], people[1]
    xt, yt, r = tree[0][0], tree[0][1], tree[1]
    x = Symbol('x')
    y = Symbol('y')
    exprs = [(x - xp) * (x - xt) + (y - yp) * (y - yt), (x - xt) ** 2 + (y - yt) ** 2 - r ** 2]
    result = solve(exprs, [x, y])
    x1, y1, x2, y2 = result[0][0], result[0][1], result[1][0], result[1][1]
    # print(x1, '\n', y1, '\n', x2, '\n', y2)
    return x1, y1, x2, y2

def coordinate(w, h, tree, people):
    #与树冠的交点坐标为(x1, y1),(x2, y2)
    x1, y1, x2, y2 = line(tree=tree, people=people)
    # print(x1, y1, x2, y2)
    xp, yp = people[0], people[1]
    xt, yt, r = tree[0][0], tree[0][1], tree[1]
    coordinates_all = []
    coordinates_select_1 = []
    # 直线1的方程为：(x - xp) / (y - yp) = (x1 - xp) / (y1 - yp)
    # 直线2的方程为：(x - xp) / (y - yp) = (x2 - xp) / (y2 - yp)
    # step1: 先求两条直线与天空的交点
    coor_x1 = (h - yp) * (x1 - xp) / (y1 - yp) + xp
    coor_x2 = (h - yp) * (x2 - xp) / (y2 - yp) + xp
    coordiante_1, coordiante_2 = [coor_x1, h], [coor_x2, h]
    coordinates_all.append(coordiante_1)
    coordinates_all.append(coordiante_2)
    # print(coordiante_1, coordiante_2)
    # step2: 求两条直线与左面墙的交点
    coor_x3 = (0 - xp) * (y1 - yp) / (x1 - xp) + yp
    coor_x4 = (0 - xp) * (y2 - yp) / (x2 - xp) + yp
    coordiante_3, coordiante_4 = [0, coor_x3], [0, coor_x4]
    coordinates_all.append(coordiante_3)
    coordinates_all.append(coordiante_4)
    # print(coordiante_3, coordiante_4)
    # step3: 求直线与右墙面的交点
    coor_x5 = (y1 - yp) * (w - xp) / (x1 - xp) + yp
    coor_x6 = (y2 - yp) * (w - xp) / (x2 - xp) + yp
    coordiante_5, coordiante_6 = [w, coor_x5], [w, coor_x6]
    coordinates_all.append(coordiante_5)
    coordinates_all.append(coordiante_6)
    # print(coordiante_5, coordiante_6)
    # print(coordinates_all)
    for i in coordinates_all:
        if 0 <= i[0] <= w and yp <= i[1] <= h:
            coordinates_select_1.append(i)
    # 对两个交点进行排序，从左往右,sort()函数先按[i][0]排然后按[i][1]排，非常符合我们的需要
    coordinates_select_1.sort()
    # print(coordinates_select_1)
    return coordinates_select_1

def view_tree(w, h, tree, people):
    '''
    Returns 返回树冠在天空方向、左墙面、右墙面方向的视角因子，即将人体与树冠的视角因子拆成三部分，用于后续辐射计算
    -------
    '''
    view_p_sky, view_p_wall1, view_p_ground, view_p_wall2 = view_wall(w=w, h=h, people=people)
    coor_1, coor_2 = coordinate(w=w, h=h, tree=tree, people=people)[0], coordinate(w=w, h=h, tree=tree, people=people)[1]
    # print(coor_1, coor_2)
    x1, y1 = 0, h
    x4, y4 = w, h
    xp, yp = people[0], people[1]
    xt, yt, r = tree[0][0], tree[0][1], tree[1]
    length_tree_p = ((xt - xp) ** 2 + (yt - yp) ** 2) ** 0.5       # 人体与树冠中心点的距离
    view_p_tree_2D = math.asin(r / length_tree_p) / np.pi             # 人与树冠的视角因子 2D
    view_p_tree = r ** 2 / (4 * ((r ** 2 + length_tree_p ** 2) ** 0.5) ** 2)   # 人与树冠的视角因子 3D
    ratio = view_p_tree / view_p_tree_2D
    # 根据光线与街道峡谷的交点得到view_p_treesky, view_p_treewall1, view_p_treewall2
    if coor_1[0] == coor_2[0] == 0:
        view_p_treesky = view_p_treewall2 = 0
        view_p_treewall1 = view_p_tree
    elif coor_1[0] == coor_2[0] == w:
        view_p_treesky = view_p_treewall1 = 0
        view_p_treewall2 = view_p_tree
    elif coor_1[1] == coor_2[1] == h:
        view_p_treewall1 = view_p_treewall2 = 0
        view_p_treesky = view_p_tree
    elif coor_1[0] == 0 and coor_2[1] == h:
        l_1 = ((xp - coor_1[0]) ** 2 + (yp - coor_1[1]) ** 2) ** 0.5
        l_2 = ((xp - x1) ** 2 + (yp - y1) ** 2) ** 0.5
        l_3 = ((xp - coor_2[0]) ** 2 + (yp - coor_2[1]) ** 2) ** 0.5
        view_p_treewall1 = math.acos((l_1 ** 2 + l_2 ** 2 - (y1 - coor_1[1]) ** 2) / (2 * l_1 * l_2)) / (2 * np.pi) * ratio
        view_p_treesky = math.acos((l_2 ** 2 + l_3 ** 2 - (coor_2[0] - x1) ** 2) / (2 * l_2 * l_3)) / (2 * np.pi) * ratio
        view_p_treewall2 = 0
    elif coor_1[1] == h and coor_2[0] == w:
        l_1 = ((xp - coor_1[0]) ** 2 + (yp - coor_1[1]) ** 2) ** 0.5
        l_2 = ((xp - x4) ** 2 + (yp - y4) ** 2) ** 0.5
        l_3 = ((xp - coor_2[0]) ** 2 + (yp - coor_2[1]) ** 2) ** 0.5
        view_p_treewall2 = math.acos((l_2 ** 2 + l_3 ** 3 - (y4 - coor_2[1]) ** 2) / (2 * l_2 * l_3)) / (2 * np.pi) * ratio
        view_p_treesky = math.acos((l_1 ** 2 + l_2 ** 2 - (x4 - coor_1[0]) ** 2) / (2 * l_1 * l_2)) / (2 * np.pi) * ratio
        view_p_treewall1 = 0
    elif coor_1[0] == 0 and coor_2[0] == w:
        l_1 = ((xp - coor_1[0]) ** 2 + (yp - coor_1[1]) ** 2) ** 0.5
        l_2 = ((xp - x1) ** 2 + (yp - y1) ** 2) ** 0.5
        l_3 = ((xp - x4) ** 2 + (yp - y4) ** 2) ** 0.5
        l_4 = ((xp - coor_2[0]) ** 2 + (yp - coor_2[1]) ** 2) ** 0.5
        view_p_treewall1 = math.acos((l_1 ** 2 + l_2 ** 2 - (y1 - coor_1[1]) ** 2) / (2 * l_1 * l_2)) / (2 * np.pi) * ratio
        view_p_treewall2 = math.acos((l_3 ** 2 + l_4 ** 2 - (y4 - coor_2[1]) ** 2) / (2 * l_3 * l_4)) / (2 * np.pi) * ratio
        view_p_treesky = view_p_sky * ratio
    print(view_p_tree, view_p_treewall1, view_p_treesky, view_p_treewall2)
    return view_p_tree, view_p_treewall1, view_p_treesky, view_p_treewall2


if __name__ == '__main__':
    a = view_tree(2, 2, tree=((1, 1.15), 0.735, 2.49, 1), people=(1.0, 0.4))


