import random
from beautifultable import BeautifulTable
import numpy as np
import time


# Гутов Віталій
# Варіант 108:
# x1_min = -30, x1_max = 0, x2_min = -35, x2_max = 10, x3_min = 0, x3_max = 20
# y_min = 200 + xc_min
# y_max = 200 + xc_max

def main():
    global m
    global y_matrix
    global average_y
    global plan_matrix_x
    global plan_matrix_for_output
    global b0, b1, b2, b3
    global xn
    global n
    x0 = 1
    x1_min = -30
    x1_max = 0
    x2_min = -35
    x2_max = 10
    x3_min = 0
    x3_max = 20

    xc_min = (x1_min + x2_min + x3_min) / 3
    xc_max = (x1_max + x2_max + x3_max) / 3

    y_min = 200 + xc_min
    y_max = 200 + xc_max

    n = 4

    xn = [[-1, -1, -1],
          [-1, 1, 1],
          [1, -1, 1],
          [1, 1, -1]]

    plan_matrix_for_output = [[x1_min, x2_min, x3_min],
                              [x1_min, x2_max, x3_max],
                              [x1_max, x2_min, x3_max],
                              [x1_max, x2_max, x3_min]]
    plan_matrix_x = [[x1_min, x1_min, x1_max, x1_max],
                     [x2_min, x2_max, x2_min, x2_max],
                     [x3_min, x3_max, x3_max, x3_min]]

    y_matrix = [[random.randint(int(y_min), int(y_max)) for _ in range(m)] for _ in range(len(xn))]

    average_y = [round(sum(i) / len(i), 3) for i in y_matrix]

    mx1 = sum(plan_matrix_x[0]) / len(plan_matrix_x[0])
    mx2 = sum(plan_matrix_x[1]) / len(plan_matrix_x[1])
    mx3 = sum(plan_matrix_x[2]) / len(plan_matrix_x[2])
    my = sum(average_y) / len(average_y)

    a1 = sum([plan_matrix_x[0][i] * average_y[i] for i in range(len(plan_matrix_x[0]))]) / len(plan_matrix_x[0])
    a2 = sum([plan_matrix_x[1][i] * average_y[i] for i in range(len(plan_matrix_x[1]))]) / len(plan_matrix_x[1])
    a3 = sum([plan_matrix_x[2][i] * average_y[i] for i in range(len(plan_matrix_x[2]))]) / len(plan_matrix_x[2])
    a11 = sum([i ** 2 for i in plan_matrix_x[0]]) / len(plan_matrix_x[0])
    a22 = sum([i ** 2 for i in plan_matrix_x[1]]) / len(plan_matrix_x[1])
    a33 = sum([i ** 2 for i in plan_matrix_x[2]]) / len(plan_matrix_x[2])
    a12 = sum([plan_matrix_x[0][i] * plan_matrix_x[1][i] for i in range(len(plan_matrix_x[0]))]) / len(plan_matrix_x[0])
    a13 = sum([plan_matrix_x[0][i] * plan_matrix_x[2][i] for i in range(len(plan_matrix_x[0]))]) / len(plan_matrix_x[0])
    a23 = sum([plan_matrix_x[1][i] * plan_matrix_x[2][i] for i in range(len(plan_matrix_x[0]))]) / len(plan_matrix_x[0])
    a32 = a23
    a31 = a13
    a21 = a12

    main_denominator = np.array([[1, mx1, mx2, mx3],
                                 [mx1, a11, a12, a13],
                                 [mx2, a12, a22, a32],
                                 [mx3, a13, a23, a33]])
    b0_numerator = np.array([[my, mx1, mx2, mx3],
                             [a1, a11, a12, a13],
                             [a2, a12, a22, a32],
                             [a3, a13, a23, a33]])
    b1_numerator = np.array([[1, my, mx2, mx3],
                             [mx1, a1, a12, a13],
                             [mx2, a2, a22, a32],
                             [mx3, a3, a23, a33]])
    b2_numerator = np.array([[1, mx1, my, mx3],
                             [mx1, a11, a1, a13],
                             [mx2, a12, a2, a32],
                             [mx3, a13, a3, a33]])

    b3_numerator = np.array([[1, mx1, mx2, my],
                             [mx1, a11, a12, a1],
                             [mx2, a12, a22, a2],
                             [mx3, a13, a23, a3]])

    b0 = np.linalg.det(b0_numerator) / np.linalg.det(main_denominator)
    b1 = np.linalg.det(b1_numerator) / np.linalg.det(main_denominator)
    b2 = np.linalg.det(b2_numerator) / np.linalg.det(main_denominator)
    b3 = np.linalg.det(b3_numerator) / np.linalg.det(main_denominator)


def kohren():
    global m
    global gp
    global s
    s = [sum([(y_matrix[j][i] - average_y[i]) ** 2 for i in range(m)]) / m for j in range(len(plan_matrix_x[0]))]
    gp = max(s) / sum(s)
    gt = 0.7679

    if gp > gt:
        m += 1
        main()


def studens():
    global n, m
    global average_y
    global b0, b1, b2, b3
    global d
    global blist
    global sb

    d = 4

    x = [[1, 1, 1, 1],
         [-1, -1, 1, 1],
         [-1, 1, -1, 1],
         [-1, 1, 1, -1]]

    sb = sum(s) / n
    s_beta_2 = sb / (n * m)
    s_beta = s_beta_2 ** (1 / 2)

    bb = [sum([average_y[i] * x[j][i] for i in range(n)]) / n for j in range(n)]
    t = [abs(bb[i]) / s_beta for i in range(n)]
    tt = 2.306
    blist = [b0, b1, b2, b3]
    for i in range(n):
        if t[i] < tt:
            blist[i] = 0
            d -= 1
    return t


def fish():
    global n
    global d
    global blist
    global plan_matrix_for_output
    global average_y
    global sb
    global fp
    y_reg = [blist[0] + blist[1] * plan_matrix_for_output[i][0] + blist[2] * plan_matrix_for_output[i][1] + blist[3] * plan_matrix_for_output[i][2] for i in range(4)]
    sad = (m / (n - d)) * int(sum([(y_reg[i] - average_y[i]) ** 2 for i in range(n)]))
    fp = sad / sb
    if fp > 4.5:
        return 'The regression equation is inadequate to the original at a significance level of 0.05'
    else:
        return 'The regression equation is adequate to the original at a significance level of 0.05'


counted_time = 0  # total program execution time
num_of_iterations = 100
for _ in range(num_of_iterations):
    m = 3

    begin_time = time.perf_counter()

    main()
    kohren()
    t_list = studens()
    fisher = fish()

    end_time = time.perf_counter()
    counted_time += end_time - begin_time

average_time = counted_time / num_of_iterations


plan_table = BeautifulTable()
headers_x = ['X{}'.format(i) for i in range(1, m+1)]
headers_y = ['Y{}'.format(i) for i in range(1, m+1)]
plan_table.columns.header = [*headers_x, *headers_y]
for i in range(len(xn)):
    plan_table.rows.append([*plan_matrix_for_output[i], *y_matrix[i]])
print('Матриця планування:')
print(plan_table)

print('Normalized equations:')
result_y = []
for i in range(len(plan_matrix_x[0])):
    result_y.append(b0 + b1 * plan_matrix_for_output[i][0] + b2 * plan_matrix_for_output[i][1] + b3 * plan_matrix_for_output[i][2])
    if round(result_y[i], 3) == round(average_y[i], 3):
        print('y = {} + {} * x1 + {} * x2 + {} * x3 = {} = {} - average y'.format(round(b0, 3), round(b1, 3), round(b2, 3), round(b3, 3), round(result_y[i], 3), average_y[i]))
    else:
        print('The obtained values do not coincide with the average values of the response function')

print('Kohren check')
print('Gp = {} < 0.7679'.format(round(gp, 3)))

print('Studens')
for i in range(len(t_list)):
    print('t{} = {}'.format(i, round(t_list[i], 3)))

print('Fisher')
print('Fp = ', round(fp, 3))
print(fisher)

print('\nAverage time of 1 iteration {} seconds'.format(average_time))
