import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
#具体指定参数
a_0=45.034/(2*np.pi)
# 求导
def f_0(y, x):
    return x**2 + y**2 - 2*x*y*np.cos(x-y) - (165)**2 / (a_0)**2

def f_1(y, x):
    return 2*y - 2*x*np.cos(x-y) - 2*x*y*np.sin(x-y)

def f_2(y, x):
    return 2 + 2*x*np.sin(x-y) + 2*x*np.sin(x-y) + 2*x*y*np.cos(x-y)


# Halley方法
def halley_method(y_initial, x, i, tol=1e-15, max_iter=100):
    y = y_initial
    for _ in range(max_iter):
        f0 = f_0(y, x)
        f1 = f_1(y, x)
        f2 = f_2(y, x)


        if i == 0:
            f0 = x**2 + y**2 - 2*x*y*np.cos(x-y) - (341-2*27.5)**2 / (a_0)**2
        
        delta_y = f0/f1*1/(1-f0*f2/(2*f1**2))
        y -= delta_y
        
        if abs(delta_y) < tol:
            break
        
        
    return y

#初始化参数
x_constant = 450/a_0
points = []
point = (a_0 * x_constant * np.cos(x_constant), a_0 * x_constant * np.sin(x_constant))
points.append(point)
y_threshold = 2 * 16 * np.pi

for i in range(0, 223):
    if i == 0:
        y_initial = x_constant + (341-2*27.5) / (a_0 * np.sqrt(x_constant**2 - 1))

    y_initial = x_constant + 159 / (a_0 * np.sqrt(x_constant**2 - 1))
    y_solution = halley_method(y_initial, x_constant, i)
    x_constant = y_solution
    
    point = (a_0 * y_solution * np.cos(y_solution), a_0 * y_solution * np.sin(y_solution))
    points.append(point)


# 可视化
plt.figure(figsize=(8, 8))

# 绘制等距螺线（最底层）
theta = np.linspace(0, 32 * np.pi, 1000)
max_radius = 55 / (2 * np.pi) * 32 * np.pi
r = np.linspace(0, max_radius, 1000)
x_helix = r * np.cos(theta)
y_helix = r * np.sin(theta)
plt.plot(x_helix, y_helix, 'lightblue', linestyle='-', linewidth=2, alpha=0.6, zorder=0)  # 设置zorder为0，使其在底层

patches = []
rectangles = []

# 绘制点
for i, point in enumerate(points):
    plt.plot(point[0], point[1], 'b', zorder=1)  

# 根据比例定义矩形
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

# 检查矩形相交
intersected_patches = []
intersection_found = False

for i in range(len(rectangles)):
    for j in range(len(rectangles)):
        if j == i or j == i-1 or j == i+1:
            continue
        path_i = Path(rectangles[i].get_xy())
        path_j = Path(rectangles[j].get_xy())
        if path_i.intersects_path(path_j):
            print(f"Intersection found between rectangles {i} and {j}")
            intersection_found = True
            intersected_patches.append((rectangles[i], i))
            intersected_patches.append((rectangles[j], j))

# 绘制矩形
p = PatchCollection(patches, facecolor='lightcoral', edgecolor='red',zorder=2)  # 设置zorder为2，使非相交矩形在点之上
plt.gca().add_collection(p)
if intersected_patches:
    p_intersected = PatchCollection([patch for patch, _ in intersected_patches], facecolor='darkred', edgecolor='red', alpha=0.8, zorder=3)  # 设置zorder为3，使相交矩形在其他矩形之上
    plt.gca().add_collection(p_intersected)

plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
