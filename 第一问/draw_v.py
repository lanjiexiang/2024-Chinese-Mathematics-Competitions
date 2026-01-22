import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

# 求导
def f_0(y, x):
    return x**2 + y**2 - 2*x*y*np.cos(x-y) - (165)**2 / (55/(2*np.pi))**2

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
            f0 = x**2 + y**2 - 2*x*y*np.cos(x-y) - (341-2*27.5)**2 / (55/(2*np.pi))**2

        delta_y = f0/f1*1/(1-f0*f2/(2*f1**2))
        y -= delta_y
        
        if abs(delta_y) < tol:
            break
        
        
    return y

# 指定秒
j_values = [0,60, 120, 180, 240, 300]
colors = ['k','b', 'g', 'r', 'c', 'm']
plt.figure(figsize=(10, 6))

for j, color in zip(j_values, colors):
    s_given = j * 100
    A = 0.5 * 55 / (2 * np.pi)
    x = 2 * np.pi * 16
    all_value = A * (x * np.sqrt(1 + x**2) + np.arcsinh(x))
    
    def f(x):
        return A * (x * np.sqrt(1 + x**2) + np.arcsinh(x)) - (all_value - s_given)

    try:
        x_solution = newton(f, 1, tol=1e-10)
    except Exception as e:
        print(f"Error for s={s_given}: {e}")
        x_solution = np.nan 

    x_constant = x_solution
    theta_coords = [x_constant]
    v_coords = [100]
    v = 100

    for i in range(223):
        y_initial = x_constant + 159 / (55 / (2 * np.pi) * np.sqrt(x_constant**2 - 1))
        y_solution = halley_method(y_initial, x_constant, i)
        x_constant = y_solution

        L = 165
        if i == 0:
            L = 341 - 2 * 27.5
        
        theta_coords.append(x_constant)
        x_1 = theta_coords[i]
        x_2 = theta_coords[i + 1]
        v = v * np.sqrt((1 + x_2**2) / (1 + x_1**2)) * (
            ((x_1 * np.sin(x_1) - np.cos(x_1)) * (x_2 * np.cos(x_2) - x_1 * np.cos(x_1)) -
             (x_1 * np.cos(x_1) + np.sin(x_1)) * (x_2 * np.sin(x_2) - x_1 * np.sin(x_1))) /
            ((x_2 * np.sin(x_2) - np.cos(x_2)) * (x_2 * np.cos(x_2) - x_1 * np.cos(x_1)) -
             (x_2 * np.cos(x_2) + np.sin(x_2)) * (x_2 * np.sin(x_2) - x_1 * np.sin(x_1)))
        )

        v_coords.append(v)

    plt.plot(range(0, len(v_coords)), v_coords, linestyle='-', color=color, label=f't={j}s')

# 调整坐标轴
plt.xlim([0, len(v_coords)])
plt.ylim([99.6, 100])

plt.xlabel('number')
plt.ylabel('v')
plt.xticks(range(0, len(v_coords), 20))
plt.yticks(np.arange(99.6, 100.02, 0.1))
plt.legend()
plt.show()
