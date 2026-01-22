import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

# 求导
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
        
        delta_y = f0/f1*1/(1-f0*f2/(2*f1**2))
        y -= delta_y
        
        if abs(delta_y) < tol:
            break
        
    return y

# 第一个点走了多远
A = 0.5 * 55 / (2 * np.pi)
x = 2 * np.pi * 16
all_value = A * (x * np.sqrt(1 + x**2) + np.arcsinh(x))
s_given = 300 * 100  # 几秒

def f(x):
    return A * (x * np.sqrt(1 + x**2) + np.log(x + np.sqrt(1 + x**2))) - (all_value - s_given)

x_solution = newton(f, 1, tol=1e-10)
x_constant = x_solution
print(x_constant)

# 迭代与可视化
points = []
point = (55 * x_constant / (2 * np.pi) * np.cos(x_constant), 55 * x_constant / (2 * np.pi) * np.sin(x_constant))
points.append(point)

for i in range(0, 223):
    if i == 0:
        y_initial = x_constant + (341 - 2 * 27.5) / (55 / (2 * np.pi) * np.sqrt(x_constant**2 - 1))

    y_initial = x_constant + 159 / (55 / (2 * np.pi) * np.sqrt(x_constant**2 - 1))
    y_solution = halley_method(y_initial, x_constant, i)
    x_constant = y_solution
    print(y_solution)
    
    point = (55 * y_solution / (2 * np.pi) * np.cos(y_solution), 55 * y_solution / (2 * np.pi) * np.sin(y_solution))
    points.append(point)

# 可视化
plt.figure(figsize=(8, 8))

# 绘制等距螺线
theta = np.linspace(0, 32 * np.pi, 1000)
max_radius = 55 / (2 * np.pi) * 32 * np.pi
r = np.linspace(0, max_radius, 1000)
x_helix = r * np.cos(theta)
y_helix = r * np.sin(theta)
plt.plot(x_helix, y_helix, 'lightblue', linestyle='-', linewidth=2, alpha=0.6)

#绘制把手位置点
for i, point in enumerate(points):
    plt.plot(point[0], point[1], 'o', color='#FFFF80') 

for i in range(1, len(points)):
    plt.plot([points[i-1][0], points[i][0]], [points[i-1][1], points[i][1]], 'r-', linewidth=2) 

# 设置坐标轴边界
plt.xlim(-1250, 1250)
plt.ylim(-1250, 1250)
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box') 
plt.show()
