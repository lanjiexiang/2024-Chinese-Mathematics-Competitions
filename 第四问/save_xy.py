import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
import pandas as pd  # 导入pandas库

# 具体指定参数
D=170
a_0 = D / (2 * np.pi)

# 求导函数定义
def f_0(y, x):
    return x**2 + y**2 - 2 * x * y * np.cos(x - y) - (165)**2 / (a_0)**2

def f_1(y, x):
    return 2 * y - 2 * x * np.cos(x - y) - 2 * x * y * np.sin(x - y)

def f_2(y, x):
    return 2 - 2*x*np.sin(x-y) - 2*x*np.sin(x-y) - 2*x*y*np.cos(x-y)

# Halley数值解
def halley_method(y_initial, x, i, tol=1e-15, max_iter=100):
    y = y_initial
    for _ in range(max_iter):
        f0 = f_0(y, x)
        f1 = f_1(y, x)
        f2 = f_2(y, x)
        if i == 0:
            f0 = x**2 + y**2 - 2 * x * y * np.cos(x - y) - (341 - 2 * 27.5)**2 / (a_0)**2
        

        delta_y = f0 / f1 * 1 / (1 - f0 * f2 / (2 * f1**2))
        y -= delta_y
        
        if abs(delta_y) < tol:
            break
        
    return y

# 两个圆心
O_1x = -2 / 3 * 170/(2*np.pi)
O_1y = 150
O_2x = 170/(2*np.pi) / 3
O_2y = -300
# 掉头角度
x_constant = 450 / a_0  
A = 0.5 * D / (2 * np.pi)
x_0 = 2 * np.pi * 16
all_value = A * (x_0 * np.sqrt(1 + x_0**2) + np.arcsinh(x_0))
#螺线长度公式
def f(x):
    return A * (x * np.sqrt(1 + x**2) + np.arcsinh(x))-(A*(x_constant * np.sqrt(1 + x_constant**2) + np.arcsinh(x_constant))+100*100)

# 起点角度与长度
x_begin = newton(f, 1, tol=1e-10)
begin_S= A*(x_constant * np.sqrt(1 + x_constant**2) + np.arcsinh(x_constant))+100*100

#两段弧长
alpha = np.arctan(450 / D*(2*np.pi))
beta = np.arctan(450 / D*(2*np.pi))
kexi_0 = np.arctan(450 / D*(2*np.pi))
R = 1 / 3 * np.sqrt(450**2 + (D/(2*np.pi))**2)
S_2 = 2 * np.pi * 2 * R * 2 * alpha / (2 * np.pi)
S_3 = 2 * np.pi * R * 2 * beta / (2 * np.pi)

#出去角度
x_out=np.pi+2*np.pi*450/D
#由距离解kexi
def judge_1(S):
    if S <= 100 * 100:
        def f(x):
            return A * (x * np.sqrt(1 + x**2) + np.arcsinh(x)) - (begin_S - S)
        x = newton(f, -20, tol=1e-10)
        kexi=-x+450*2*np.pi/D
    elif S >= 100 * 100 and S <= 100*100 +S_2:
        alpha_0= (S - 100 * 100)/(2*R) 
        kexi=2*alpha_0/3
    elif S >= 100*100 +S_2 and S <= 100*100 + S_2+S_3:
        beta_0=(S - 100 * 100-S_2)/(R)
        kexi=(beta_0+4*kexi_0)/3
    else:
        
        S_out=A * ((x_out-np.pi) * np.sqrt(1 + (x_out-np.pi)**2) + np.arcsinh(x_out-np.pi))
        def f(x):
            return A * ((x-np.pi) * np.sqrt(1 + (x-np.pi)**2) + np.arcsinh(x-np.pi))-S_out-(S-S_2-S_3-100*100)
        x= newton(f, 1, tol=1e-10)
        kexi=x+2*kexi_0-np.pi-(2*np.pi*450)/D
    return kexi
#由kexi解位置
def distance(kexi):
    if kexi <= 0 :
        x = -kexi + 450/(D/(2*np.pi))
        x_coord= a_0 * x * np.cos(x)
        y_coord= a_0 * x * np.sin(x)
    elif kexi <= 4/3*kexi_0:
        theta = x_constant -2*2*np.pi
        rotation_angle = theta - 0.5*np.pi
        alpha_0=3/2*kexi
        delta = alpha - alpha_0
        point = (2 * R * np.cos(delta) + O_1x, 2 * R * np.sin(delta) + O_1y)
        rotation_matrix = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle)],
                            [np.sin(rotation_angle),  np.cos(rotation_angle)]])
        (x_coord,y_coord) = np.dot(rotation_matrix, point)
    elif kexi <= 2*kexi_0:
        theta = x_constant -2*2*np.pi
        rotation_angle = theta - 0.5*np.pi
        beta_0=3*kexi-4*kexi_0
        delta = beta - beta_0
        point = (-R * np.cos(delta) + O_2x, R * np.sin(delta) + O_2y)
        rotation_matrix = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle)],
                            [np.sin(rotation_angle),  np.cos(rotation_angle)]])
        (x_coord,y_coord) = np.dot(rotation_matrix, point)
    else :
        x=kexi-2*kexi_0+np.pi+2*np.pi*450/D
        x_coord=D/(2*np.pi)*(x-np.pi)*np.cos(x)
        y_coord=D/(2*np.pi)*(x-np.pi)*np.sin(x)

    return (x_coord, y_coord)


# 初始化变量
all_x = []
all_y = [] 
#遍历时间-100s认为是0s
for t in np.arange(0, 201):
    S = 100 * t
    kexi = judge_1(S)
    point = distance(kexi)
    points = [(point, kexi)]
    kexi_old = kexi

    for i in np.arange(0, 223):
        if i == 0:
            L = 341 - 55
        else:
            L = 165
        (x_old, y_old) = point

        def f(xi):
            x_new, y_new = distance(xi)
            return (x_new - x_old)**2 + (y_new - y_old)**2 - L**2

        try:
            # 自适应调整初值
            if kexi_old > 2 * kexi_0:
                solve_kexi = newton(f, kexi_old - 0.5, tol=1e-3, maxiter=500)

                # 检查解是否比旧解小且不同
                if -solve_kexi + kexi_old < np.pi and abs(solve_kexi - kexi_old) > 1e-6:
                    point = distance(solve_kexi)
                    kexi_old = solve_kexi
                    points.append((point, solve_kexi))
            else:
                solve_kexi = newton(f, kexi_old - 1.2, tol=1e-3, maxiter=500)

                # 检查解是否比旧解小且不同
                if -solve_kexi + kexi_old < np.pi and abs(solve_kexi - kexi_old) > 1e-6:
                    point = distance(solve_kexi)
                    kexi_old = solve_kexi
                    points.append((point, solve_kexi))
        except RuntimeError:
            print("no")

    # 从 points 中提取所有的 x 和 y 坐标
    x_coords = [p[0][0] for p in points]
    y_coords = [p[0][1] for p in points]
    all_x.append(x_coords)
    all_y.append(y_coords)

# 将 x 和 y 坐标分别保存到两个 DataFrame 中，并写入 Excel 文件
df_x = pd.DataFrame(all_x).transpose()  
df_y = pd.DataFrame(all_y).transpose() 
with pd.ExcelWriter('points_coordinates.xlsx') as writer:
    df_x.to_excel(writer, sheet_name='X_Coordinates', index=False)
    df_y.to_excel(writer, sheet_name='Y_Coordinates', index=False)

print("x 和 y 坐标已保存到 points_coordinates.xlsx 文件中")

# 绘制所有点
for point, kexi in points:
    if kexi < 0:
        plt.plot(point[0], point[1], 'ko')  # 黑色点
    elif 0 <= kexi < 4 / 3 * kexi_0:
        plt.plot(point[0], point[1], 'bo')  # 蓝色点
    elif 4 / 3 * kexi_0 <= kexi < 2 * kexi_0:
        plt.plot(point[0], point[1], 'go')  # 绿色点
    else:
        plt.plot(point[0], point[1], 'yo')  # 黄色点

# 设置图形属性
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()