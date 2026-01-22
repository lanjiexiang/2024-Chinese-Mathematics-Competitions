import numpy as np
from scipy.optimize import newton
from matplotlib.patches import Polygon
from matplotlib.path import Path

#求导
def f_0(y, x):
    return x**2 + y**2 - 2*x*y*np.cos(x-y) - (165)**2 / (55/(2*np.pi))**2

def f_1(y, x):
    return 2*y - 2*x*np.cos(x-y) - 2*x*y*np.sin(x-y)

def f_2(y, x):
    return 2 - 2*x*np.sin(x-y) - 2*x*np.sin(x-y) - 2*x*y*np.cos(x-y)


# Halley方法
def halley_method(y_initial, x, i, tol=1e-15, max_iter=100):
    y = y_initial
    for _ in range(max_iter):
        f0 = f_0(y, x)
        f1 = f_1(y, x)
        f2 = f_2(y, x)

        if i == 0:
            f0 = x**2 + y**2 - 2*x*y*np.cos(x-y) - (341-2*27.5)**2 / (55/(2*np.pi))**2
        
        delta_y =f0/f1*1/(1-f0*f2/(2*f1**2))
        y -= delta_y
        
        if abs(delta_y) < tol:
            break
        
    return y

# 初始化参数
A = 0.5 * 55 / (2 * np.pi)
x = 2 * np.pi * 16
all = A * (x * np.sqrt(1 + x**2) + np.log(x + np.sqrt(1 + x**2)))

#开始检查，j为时间
for j in np.arange(100, 500, 0.01):
    s_given = j * 100
    print(j)

    def f(x):
        return A * (x * np.sqrt(1 + x**2) + np.arcsinh(x)) - (all - s_given)

    x_solution = newton(f, 1, tol=1e-10)
    x_constant = x_solution


    points = []
    point = (55 * x_constant / (2 * np.pi) * np.cos(x_constant), 55 * x_constant / (2 * np.pi) * np.sin(x_constant))
    points.append(point)
    y_threshold = 2 * 16 * np.pi

    for i in range(0, 223):
        if i == 0:
            y_initial = x_constant + (341-2*27.5) / (55 / (2 * np.pi) * np.sqrt(x_constant**2 - 1))

        y_initial = x_constant + 159 / (55 / (2 * np.pi) * np.sqrt(x_constant**2 - 1))
        y_solution = halley_method(y_initial, x_constant, i)

        # 是否进入螺旋
        if np.isclose(y_solution, y_threshold, atol=1e-10):
            print(f"Stopping iteration at i={i} because y_solution reached the threshold {y_threshold}")
            break

        x_constant = y_solution

        point = (55 * y_solution / (2 * np.pi) * np.cos(y_solution), 55 * y_solution / (2 * np.pi) * np.sin(y_solution))
        points.append(point)


    patches = []
    rectangles = []
    #根据比例尺画矩形
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


    intersected_patches = []
    intersection_found = False 
    #判断矩形相交
    for i in range(len(rectangles)):
        for j_rect in range(len(rectangles)):
            if j_rect == i or j_rect == i-1 or j_rect == i+1: 
                continue
            path_i = Path(rectangles[i].get_xy())
            path_j = Path(rectangles[j_rect].get_xy())
            if path_i.intersects_path(path_j):
                print(f"Collision detected at j = {j}, between rectangles {i} and {j_rect}")
                intersection_found = True
                intersected_patches.append(rectangles[i])
                intersected_patches.append(rectangles[j_rect])

                break
        if intersection_found:
            break

    # 一旦矩形相交结束
    if intersection_found:
        break 