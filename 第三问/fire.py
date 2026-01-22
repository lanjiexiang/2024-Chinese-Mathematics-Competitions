import numpy as np
from matplotlib.patches import Polygon
from matplotlib.path import Path
import matplotlib.pyplot as plt


# 求导
def f_0(y, x, a):
    return x**2 + y**2 - 2*x*y*np.cos(x-y) - (165)**2 / (a)**2

def f_1(y, x):
    return 2*y - 2*x*np.cos(x-y) - 2*x*y*np.sin(x-y)

def f_2(y, x):
    return 2 - 2*x*np.sin(x-y) - 2*x*np.sin(x-y) + 2*x*y*np.cos(x-y)


# Halley方法
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

# 定义优化目标函数
def check(k,j):
    print(k,j)
    no_intersection_found = True  
    a_0 = k / (2 * np.pi)
    x_constant = j
    points = []
    point = (a_0 * x_constant * np.cos(x_constant), a_0 * x_constant * np.sin(x_constant))
    points.append(point)
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
            for j in range(len(rectangles)):
                if j == i or j == i-1 or j == i+1:
                    continue
                path_i = Path(rectangles[i].get_xy())
                path_j = Path(rectangles[j].get_xy())
                if path_i.intersects_path(path_j):
                    intersection_found = True
                    no_intersection_found = False  
                    break
            if intersection_found:
                break


    if no_intersection_found:
        return np.sqrt(k)
    else:
        return 1000000


# 模拟退火参数
initial_temperature = 1000
final_temperature = 1
alpha = 0.95  # 温度衰减率
max_iterations = 300  # 最大迭代轮数
k_bounds = (45, 46)

# 初始化种群
def initialize_individual():
    k = np.random.uniform(*k_bounds)
    j_lower_bound = 450 * 2 * np.pi / k
    j_upper_bound = 450 * 2 * np.pi / k + 3 * np.pi
    j = np.random.uniform(j_lower_bound, j_upper_bound)
    return k, j

# 评估适应度
def evaluate_individual(individual):
    k, j = individual
    return check(k, j)

# 生成邻域解
def generate_neighbor(current_solution):
    k_new = np.random.uniform(*k_bounds)
    j_lower_bound = 450 * 2 * np.pi / k_new
    j_upper_bound = 450 * 2 * np.pi / k_new + 3 * np.pi
    j_new = np.random.uniform(j_lower_bound, j_upper_bound)
    return k_new, j_new

# 模拟退火算法
def simulated_annealing(initial_temp, final_temp, alpha, max_iter):
    current_solution = initialize_individual()
    current_fitness = evaluate_individual(current_solution)
    best_solution = current_solution
    best_fitness = current_fitness
    temperature = initial_temp
    fitness_history = []

    for _ in range(max_iter):
        neighbor = generate_neighbor(current_solution)
        neighbor_fitness = evaluate_individual(neighbor)
        delta_fitness = neighbor_fitness - current_fitness

        if delta_fitness < 0 or np.random.rand() < np.exp(-delta_fitness / temperature):
            current_solution = neighbor
            current_fitness = neighbor_fitness

            if current_fitness < best_fitness:
                best_solution = current_solution
                best_fitness = current_fitness

        temperature *= alpha
        fitness_history.append(best_fitness)

        if temperature < final_temp:
            break

    return best_solution, best_fitness, fitness_history

# 运行模拟退火算法
best_solution, best_value, fitness_history = simulated_annealing(
    initial_temperature, final_temperature, alpha, max_iterations
)

print("最佳解决方案:", best_solution)
print("最佳函数值:", best_value)
k_best, j_best = best_solution
print("最佳解决方案 (k, j):", (k_best, j_best))
# 可视化学习曲线
plt.figure(figsize=(8, 6))
plt.plot(range(len(fitness_history)), fitness_history, 'b-', label='Best Fitness')
plt.title('Simulated Annealing Learning Curve')
plt.xlabel('Iteration')
plt.ylabel('Best Fitness')
plt.legend()
plt.show()