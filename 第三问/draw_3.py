import numpy as np
import matplotlib.pyplot as plt
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

# 存储 k 和 j 的值以及对应的 flag
results = []


for k in np.arange(45.0329,45.0340,0.0001):
        print(f"Testing k = {k}")
        
        for j in np.arange(450 / (k / (2 * np.pi))+0.96, 450 / (k / (2 * np.pi))+1.05,0.001):
            print(f"Testing j = {j}")
            
            no_intersection_found = True  # This will track if all iterations for this k have no intersections
            
            a_0 = k / (2 * np.pi)
            x_constant = j
            points = []
            point = (a_0 * x_constant * np.cos(x_constant), a_0 * x_constant * np.sin(x_constant))
            points.append(point)
            y_threshold = 2 * 16 * np.pi

            for i in range(0, 223):
                if i == 0:
                    y_initial = x_constant + (341 - 2 * 27.5) / (a_0 * np.sqrt(x_constant**2 - 1))

                y_initial = x_constant + 159 / (a_0 * np.sqrt(x_constant**2 - 1))
                y_solution = halley_method(y_initial, x_constant, i, a_0)
                x_constant = y_solution
                point = (a_0 * y_solution * np.cos(y_solution), a_0 * y_solution * np.sin(y_solution))
                points.append(point)

            patches = []
            rectangles = []



            # 根据坐标比例画矩形
            for i in range(1, len(points)):
                x1, y1 = points[i-1]
                x2, y2 = points[i]
                
                length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                width = length * 30 / (341-55) if i == 1 else length * 30 / 165
                extension = length * 27.5 / (341-55) if i == 1 else length * 27.5 / 165

                dx = x2 - x1
                dy = y2 - y1
                norm = np.sqrt(dx**2 + dy**2)
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

            # 检查相交，相邻不检查
            intersection_found = False

            for i in range(len(rectangles)):
                for t in range(len(rectangles)):
                    if t == i or t == i-1 or t == i+1:
                        continue
                    path_i = Path(rectangles[i].get_xy())
                    path_t = Path(rectangles[t].get_xy())
                    if path_i.intersects_path(path_t):
                        print(f"Intersection found between rectangles {i} and {t}")
                        intersection_found = True
                        no_intersection_found = False  # Set this to False since we found an intersection
                        break
                if intersection_found:
                    break

            # 设置 flag
            flag = 0 if intersection_found else 1
            results.append((k, j-450 / (k / (2 * np.pi)), flag))

# 绘制 k 和 j 的图
plt.figure(figsize=(8, 6))
for k, j, flag in results:
    if flag == 0:
        plt.plot(j, k, 'ro')  # 红色点表示相交
    else:
        plt.plot(j, k, 'bo')  # 蓝色点表示不相交

plt.xlabel('j')
plt.ylabel('k')
plt.grid()
plt.show()