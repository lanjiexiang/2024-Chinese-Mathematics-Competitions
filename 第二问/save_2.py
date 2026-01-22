import numpy as np
from scipy.optimize import newton
import pandas as pd
import math

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


        # 第一个板凳长度不一样
        if i == 0:
            f0 = x**2 + y**2 - 2*x*y*np.cos(x-y) - (341-2*27.5)**2 / (55/(2*np.pi))**2

        delta_y = f0/f1*1/(1-f0*f2/(2*f1**2))
        y -= delta_y
        
        if abs(delta_y) < tol:
            break
        
    return y


# 创建两个空的 DataFrame 用于存储数据
x_results_df = pd.DataFrame()
y_results_df = pd.DataFrame()
v_results_df = pd.DataFrame()

# 遍历 s 的不同值
s_given = 412.47383 * 100

# 计算 x 的初始值
A = 0.5 * 55 / (2 * np.pi)
x = 2 * np.pi * 16
all_value = A * (x * np.sqrt(1 + x**2) + np.arcsinh(x))
    
def f(x):
    return A * (x * np.sqrt(1 + x**2) + np.log(x + np.sqrt(1 + x**2))) - (all_value - s_given)

try:
    x_solution = newton(f, 1, tol=1e-10)
except Exception as e:
    print(f"Error for s={s_given}: {e}")
    x_solution = np.nan 

x_constant = x_solution
x_coords = [55 * x_constant / (2 * np.pi) * np.cos(x_constant)]  
y_coords = [55 * x_constant / (2 * np.pi) * np.sin(x_constant)]  
theta_coords = [x_constant]
v_coords=[100]
v=100

for i in range(223):
    y_initial = x_constant + 159 / (55 / (2 * np.pi) * np.sqrt(x_constant**2 - 1))
    y_solution = halley_method(y_initial, x_constant, i)
    x_constant = y_solution
    L=165
    if i == 0 :
        L=341-2*27.5
    a=55/2/np.pi        
    # 记录坐标
    theta_coords.append(x_constant)
    x_1=theta_coords[i]
    x_2=theta_coords[i+1]
    v = v * np.sqrt((1 + x_2**2) / (1 + x_1**2)) * (
    ((x_1 * np.sin(x_1) - np.cos(x_1)) * (x_2 * np.cos(x_2) - x_1 * np.cos(x_1)) -
    (x_1 * np.cos(x_1) + np.sin(x_1)) * (x_2 * np.sin(x_2) - x_1 * np.sin(x_1))) /
    ((x_2 * np.sin(x_2) - np.cos(x_2)) * (x_2 * np.cos(x_2) - x_1 * np.cos(x_1)) -
    (x_2 * np.cos(x_2) + np.sin(x_2)) * (x_2 * np.sin(x_2) - x_1 * np.sin(x_1))))

    v_coords.append(v)
    # 记录坐标
    x_coord = 55 * y_solution / (2 * np.pi) * np.cos(y_solution)
    y_coord = 55 * y_solution / (2 * np.pi) * np.sin(y_solution)
    x_coords.append(x_coord)
    y_coords.append(y_coord)
        
# 将结果保存到 DataFrame
x_results_df[f'{math.ceil((s_given-1)/100)}s'] = x_coords
y_results_df[f'{math.ceil((s_given-1)/100)}s'] = y_coords
v_results_df[f'{math.ceil((s_given-1)/100)}s'] = v_coords

# 保存到 Excel 文件
x_results_df.to_excel('x_2.xlsx', index=True)
y_results_df.to_excel('y_2.xlsx', index=True)
v_results_df.to_excel('v_2.xlsx', index=True)

