import numpy as np
from scipy.optimize import newton
import pandas as pd  

# 具体指定参数
D = 170
a_0 = D / (2 * np.pi)

# 求导函数定义
def f_0(y, x):
    return x**2 + y**2 - 2 * x * y * np.cos(x - y) - (165)**2 / (a_0)**2

def f_1(y, x):
    return 2 * y - 2 * x * np.cos(x - y) - 2 * x * y * np.sin(x - y)

def f_2(y, x):
    return 2 - 2 * x * np.sin(x - y) - 2 * x * np.sin(x - y) + 2 * x * y * np.cos(x - y)


# Halley方法
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
# 螺线长度公式
def f(x):
    return A * (x * np.sqrt(1 + x**2) + np.arcsinh(x))-(A*(x_constant * np.sqrt(1 + x_constant**2) + np.arcsinh(x_constant))+100*100)

# 起点角度与长度
x_begin = newton(f, 1, tol=1e-10)
begin_S = A * (x_constant * np.sqrt(1 + x_constant**2) + np.arcsinh(x_constant)) + 100*100

# 两段弧长
alpha = np.arctan(450 / D * (2*np.pi))
beta = np.arctan(450 / D * (2*np.pi))
kexi_0 = np.arctan(450 / D * (2*np.pi))
R = 1 / 3 * np.sqrt(450**2 + (D/(2*np.pi))**2)
S_2 = 2 * np.pi * 2 * R * 2 * alpha / (2 * np.pi)
S_3 = 2 * np.pi * R * 2 * beta / (2 * np.pi)

# 出去角度
x_out = np.pi + 2*np.pi*450/D

#由S解kexi
def judge_1(S):
    if S <= 100 * 100:
        def f(x):
            return A * (x * np.sqrt(1 + x**2) + np.arcsinh(x)) - (begin_S - S)
        x = newton(f, -20, tol=1e-10)
        kexi = -x + 450*2*np.pi/D
    elif S >= 100 * 100 and S <= 100*100 + S_2:
        alpha_0 = (S - 100 * 100)/(2*R) 
        kexi = 2*alpha_0/3
    elif S >= 100*100 + S_2 and S <= 100*100 + S_2 + S_3:
        beta_0 = (S - 100 * 100 - S_2)/(R)
        kexi = (beta_0 + 4*kexi_0)/3
    else:
        S_out = A * ((x_out-np.pi) * np.sqrt(1 + (x_out-np.pi)**2) + np.arcsinh(x_out-np.pi))
        def f(x):
            return A * ((x-np.pi) * np.sqrt(1 + (x-np.pi)**2) + np.arcsinh(x-np.pi)) - S_out - (S-S_2-S_3-100*100)
        x = newton(f, 1, tol=1e-10)
        kexi = x + 2*kexi_0 - np.pi - (2*np.pi*450)/D
    return kexi

#由kexi解位置
def distance(kexi):
    if kexi <= 0:
        x = -kexi + 450/(D/(2*np.pi))
        x_coord = a_0 * x * np.cos(x)
        y_coord = a_0 * x * np.sin(x)
    elif kexi <= 4/3*kexi_0:
        theta = x_constant - 2*2*np.pi
        rotation_angle = theta - 0.5*np.pi
        alpha_0 = 3/2*kexi
        delta = alpha - alpha_0
        point = (2 * R * np.cos(delta) + O_1x, 2 * R * np.sin(delta) + O_1y)
        rotation_matrix = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle)],
                            [np.sin(rotation_angle),  np.cos(rotation_angle)]])
        (x_coord, y_coord) = np.dot(rotation_matrix, point)
    elif kexi <= 2*kexi_0:
        theta = x_constant - 2*2*np.pi
        rotation_angle = theta - 0.5*np.pi
        beta_0 = 3*kexi - 4*kexi_0
        delta = beta - beta_0
        point = (-R * np.cos(delta) + O_2x, R * np.sin(delta) + O_2y)
        rotation_matrix = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle)],
                            [np.sin(rotation_angle),  np.cos(rotation_angle)]])
        (x_coord, y_coord) = np.dot(rotation_matrix, point)
    else:
        x = kexi - 2*kexi_0 + np.pi + 2*np.pi*450/D
        x_coord = D/(2*np.pi)*(x-np.pi)*np.cos(x)
        y_coord = D/(2*np.pi)*(x-np.pi)*np.sin(x)

    return (x_coord, y_coord)
#由kexi解速度方向
def direction(kexi):
    if kexi <= 0:
        x = -kexi + 450/(D/(2*np.pi))
        v_x = 1/np.sqrt(x**2+1)*(x*np.sin(x)-np.cos(x))
        v_y = 1/np.sqrt(x**2+1)*(-x*np.cos(x)-np.sin(x))
    elif kexi <= 4/3*kexi_0:
        alpha_0 = 3/2*kexi
        v_x = np.cos(2*np.pi*450/D-np.pi+kexi_0-alpha_0)
        v_y = np.sin(2*np.pi*450/D-np.pi+kexi_0-alpha_0)
    elif kexi <= 2*kexi_0:
        beta_0 = 3*kexi - 4*kexi_0
        v_x = np.cos(2*np.pi*450/D-np.pi-kexi_0+beta_0)
        v_y = np.sin(2*np.pi*450/D-np.pi-kexi_0+beta_0)
    else:
        x = kexi - 2*kexi_0 + np.pi + 2*np.pi*450/D
        v_x = 1/np.sqrt(1+(x-np.pi)**2)*((x-np.pi)*np.sin(x)-np.cos(x))
        v_y = 1/np.sqrt(1+(x-np.pi)**2)*(-(x-np.pi)*np.cos(x)-np.sin(x))
    
    return (-v_x, -v_y)

# 存储速度数据
all_Vs = []

for t in np.arange(0, 201):
    print(t)
    S = 100 * t
    kexi = judge_1(S)
    point = distance(kexi)
    points = [(point, kexi)]
    kexi_old = kexi
    V = 100
    (v_xold, v_yold) = direction(kexi)
    V_xs = [v_xold * V / np.sqrt(v_xold**2 + v_yold**2)]
    V_ys = [v_yold * V / np.sqrt(v_xold**2 + v_yold**2)]
    V_s=[np.sqrt((v_xold * V / np.sqrt(v_xold**2 + v_yold**2))**2+(v_yold * V / np.sqrt(v_xold**2 + v_yold**2))**2)]
    
    for i in np.arange(0, 223):
        if i == 0:
            L = 341 - 55
        else:
            L = 165
        
        # 上一个点
        (x_old, y_old) = point

        def f(xi):
            x_new, y_new = distance(xi)
            return (x_new - x_old)**2 + (y_new - y_old)**2 - L**2

        try:#自适应调整迭代初值
            if kexi_old > 2 * kexi_0:
                solve_kexi = newton(f, kexi_old - 0.5, tol=1e-3, maxiter=500)
                if -solve_kexi + kexi_old < np.pi and abs(solve_kexi - kexi_old) > 1e-6:
                    point = distance(solve_kexi)
                    kexi_old = solve_kexi
                    points.append((point, solve_kexi))
            else:
                solve_kexi = newton(f, kexi_old - 1.2, tol=1e-3, maxiter=500)
                if -solve_kexi + kexi_old < np.pi and abs(solve_kexi - kexi_old) > 1e-6:
                    point = distance(solve_kexi)
                    kexi_old = solve_kexi
                    points.append((point, solve_kexi))

            # 计算速度
            (v_x, v_y) = direction(solve_kexi)
            # 上一个点和现在点
            rx_old = x_old
            ry_old = y_old
            (rx_new, ry_new) = point
            l_x = rx_old - rx_new
            l_y = ry_old - ry_new
            # 速度大小递推
            V = V * (v_xold * l_x + v_yold * l_y) / (v_x * l_x + v_y * l_y)
            v_xold = v_x
            v_yold = v_y
            V_xs.append(v_x * V / np.sqrt(v_x**2 + v_y**2))
            V_ys.append(v_y * V / np.sqrt(v_x**2 + v_y**2))
            V_s.append(abs(V))

        except RuntimeError:
            print("no")


    all_Vs.append(V_s)


# 将速度数据保存到Excel文件
df_x = pd.DataFrame(all_Vs).T  
df_x.to_excel('V.xlsx', index=True, header=False)


