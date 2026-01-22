import numpy as np
from matplotlib.patches import Polygon
from matplotlib.path import Path

# 求导
def f_0(y, x, a):
    return x**2 + y**2 - 2*x*y*np.cos(x-y) - (165)**2 / (a)**2

def f_1(y, x):
    return 2*y - 2*x*np.cos(x-y) - 2*x*y*np.sin(x-y)

def f_2(y, x):
    return 2 + 2*x*np.sin(x-y) + 2*x*np.sin(x-y) + 2*x*y*np.cos(x-y)


# Halley数值解
def halley_method(y_initial, x, i, a, tol=1e-15, max_iter=100):
    y = y_initial
    for _ in range(max_iter):
        f0 = f_0(y, x, a)
        f1 = f_1(y, x)
        f2 = f_2(y, x)


        if i == 0:
            f0 = x**2 + y**2 - 2*x*y*np.cos(x-y) - (341-2*27.5)**2 / (a)**2
        

        delta_y = f0/f1*1/(1-f0*f2/(2*f1**2))
        y -= delta_y
        
        if abs(delta_y) < tol:
            break
        
    return y

def adaptive_grid_search():
    # 由模拟退火算法筛选数据后得到的初始查询范围
    initial_k_range = (45.03, 45.04)
    d1 = 0.001
    k_step = d1
    tol_k = 1e-4
    k_values = np.arange(initial_k_range[0], initial_k_range[1], k_step)
    count = 1

    # 当螺距的精度达到6位小数时停止划分网格
    while d1 >= tol_k:
        results = []
        print("第", count, "轮自适应网格")
        count += 1

        for k in k_values:
            print(f"Testing k = {k}")
            
            # 使用当前的 k 值初始化 j_values
            j_values = np.arange(450 / (k / (2 * np.pi)) + 0.96, 450 / (k / (2 * np.pi)) + 1.01, 0.001)

            for j in j_values:
                print(f"Testing j = {j}")
                # 开始查询，判断给定 k, j 下是否发生碰撞

                a_0 = k / (2 * np.pi)
                x_constant = j
                points = []
                point = (a_0 * x_constant * np.cos(x_constant), a_0 * x_constant * np.sin(x_constant))
                points.append(point)
                
                # 求解此时刻队伍位置
                for i in range(0, 223):
                    if i == 0:
                        y_initial = x_constant + (341 - 2 * 27.5) / (a_0 * np.sqrt(x_constant ** 2 - 1))

                    y_initial = x_constant + 159 / (a_0 * np.sqrt(x_constant ** 2 - 1))
                    y_solution = halley_method(y_initial, x_constant, i, a_0)
                    x_constant = y_solution
                    point = (a_0 * y_solution * np.cos(y_solution), a_0 * y_solution * np.sin(y_solution))
                    points.append(point)

                patches = []
                rectangles = []
                # 根据坐标比例画矩形
                for i in range(1, len(points)):
                    x1, y1 = points[i - 1]
                    x2, y2 = points[i]

                    length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
                    width = length * 30 / (341 - 55) if i == 1 else length * 30 / 165
                    extension = length * 27.5 / (341 - 55) if i == 1 else length * 27.5 / 165

                    dx = x2 - x1
                    dy = y2 - y1
                    norm = np.sqrt(dx ** 2 + dy ** 2)
                    dx /= norm
                    dy /= norm

                    perp_dx = -dy
                    perp_dy = dx

                    corner1 = (x1 - dx * extension + perp_dx * width / 2, y1 - dy * extension + perp_dy * width / 2)
                    corner2 = (x1 - dx * extension - perp_dx * width / 2, y1 - dy * extension - perp_dy * width / 2)
                    corner3 = (x2 + dx * extension + perp_dx * width / 2, y2 + dy * extension + perp_dy * width / 2)
                    corner4 = (x2 + dx * extension - perp_dx * width / 2, y2 + dy * extension - perp_dy * width / 2)

                    rectangle = Polygon([corner1, corner2, corner4, corner3], closed=True)
                    patches.append(rectangle)
                    rectangles.append(rectangle)

                # 矩形相交检验
                intersection_found = False
                for i in range(len(rectangles)):
                    for t in range(len(rectangles)):
                        if t == i or t == i - 1 or t == i + 1:
                            continue
                        path_i = Path(rectangles[i].get_xy())
                        path_t = Path(rectangles[t].get_xy())
                        if path_i.intersects_path(path_t):
                            print(f"Intersection found between rectangles {i} and {t}")
                            intersection_found = True
                            break
                    if intersection_found:
                        break

                flag = 0 if intersection_found else 1
                results.append((k, j - 450 / (k / (2 * np.pi)), flag))

        # 检验完毕后，找到红色点中螺距最大的点，更新新一轮的网格
        max_k_with_intersection = max((k for k, j, flag in results if flag == 0), default=None)
        # 螺距上下界更新为上下变换一格，十等分
        red_points = [(k, j) for k, j, flag in results if flag == 0]
        max_red_point_j = max(red_points, key=lambda x: x[0])[1]
        d1 /= 2
        # 时间上界更新为右侧第一个蓝点的时间，下界更新为左侧第一个蓝点的时间，50等分
        left_limit = next((j for k, j, flag in results if flag == 1 and j < max_red_point_j), max_red_point_j - 0.05)
        right_limit = next((j for k, j, flag in results if flag == 1 and j > max_red_point_j), max_red_point_j + 0.05)
        j_step = (right_limit - left_limit) / 50
        j_values = np.arange(left_limit, right_limit, j_step)
        k_values = np.arange(max_k_with_intersection - d1, max_k_with_intersection + d1, d1)

    print("最大螺距为：", max_k_with_intersection)


adaptive_grid_search()
