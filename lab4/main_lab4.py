import random
from beautifultable import BeautifulTable
from scipy.stats import f, t

# Гутов Віталій
# Варіант 108:
# x1_min = -5, x1_max = 15,
# x2_min = -15, x2_max = 35,
# x3_min = 15, x3_max = 30
# y_min = 200 + xc_min
# y_max = 200 + xc_max


def main():
    global x1_min, x1_max, x2_min, x2_max, x3_min, x3_max
    global m
    global y_matrix
    global average_y
    global n
    global b0, b1, b2, b3, b12, b13, b23, b123
    global plan_matrix, plan_matrix_normal
    global f1, f2, f3, q1, q
    m = 3
    n = 8

    f1 = m-1
    f2 = n
    f3 = f1 * f2
    q = 0.05
    q1 = q / f2

    x1_min = -5
    x1_max = 15
    x2_min = -15
    x2_max = 35
    x3_min = 15
    x3_max = 30

    xc_min = (x1_min + x2_min + x3_min) / 3
    xc_max = (x1_max + x2_max + x3_max) / 3

    y_min = 200 + xc_min
    y_max = 200 + xc_max

    y_matrix = [[random.randint(int(y_min), int(y_max)) for _ in range(m)] for _ in range(n)]

    average_y = [round(sum(i) / len(i), 3) for i in y_matrix]

    norm_x = [[-1, -1, -1],
              [-1, -1, 1],
              [-1, 1, -1],
              [-1, 1, 1],
              [1, -1, -1],
              [1, -1, 1],
              [1, 1, -1],
              [1, 1, 1]]

    b0 = sum(average_y) / n
    b1 = sum([average_y[i] * norm_x[i][0] for i in range(n)]) / n
    b2 = sum([average_y[i] * norm_x[i][1] for i in range(n)]) / n
    b3 = sum([average_y[i] * norm_x[i][2] for i in range(n)]) / n
    b12 = sum([average_y[i] * norm_x[i][0] * norm_x[i][1] for i in range(n)]) / n
    b13 = sum([average_y[i] * norm_x[i][0] * norm_x[i][2] for i in range(n)]) / n
    b23 = sum([average_y[i] * norm_x[i][1] * norm_x[i][2] for i in range(n)]) / n
    b123 = sum([average_y[i] * norm_x[i][0] * norm_x[i][1] * norm_x[i][2] for i in range(n)]) / n

    plan_matrix = [[x1_min, x2_min, x3_min, x1_min * x2_min, x1_min * x3_min, x2_min * x3_min, x1_min * x2_min * x3_min],
                   [x1_min, x2_min, x3_max, x1_min * x2_min, x1_min * x3_max, x2_min * x3_max, x1_min * x2_min * x3_max],
                   [x1_min, x2_max, x3_min, x1_min * x2_max, x1_min * x3_min, x2_max * x3_min, x1_min * x2_max * x3_min],
                   [x1_min, x2_max, x3_max, x1_min * x2_max, x1_min * x3_max, x2_max * x3_max, x1_min * x2_max * x3_max],
                   [x1_max, x2_min, x3_min, x1_max * x2_min, x1_max * x3_min, x2_min * x3_min, x1_max * x2_min * x3_min],
                   [x1_max, x2_min, x3_max, x1_max * x2_min, x1_max * x3_max, x2_min * x3_max, x1_max * x2_min * x3_max],
                   [x1_max, x2_max, x3_min, x1_max * x2_max, x1_max * x3_min, x2_max * x3_min, x1_max * x2_max * x3_min],
                   [x1_max, x2_max, x3_max, x1_max * x2_max, x1_max * x3_max, x2_max * x3_max, x1_max * x2_max * x3_max]]

    plan_matrix_normal = [[-1, -1, -1, 1, 1, 1, -1],
                          [-1, -1, 1, 1, -1, -1, 1],
                          [-1, 1, -1, -1, 1, -1, 1],
                          [-1, 1, 1, -1, -1, 1, -1],
                          [1, -1, -1, -1, -1, 1, 1],
                          [1, -1, 1, -1, 1, -1, -1],
                          [1, 1, -1, 1, -1, -1, -1],
                          [1, 1, 1, 1, 1, 1, 1]]

    result_y = []
    for i in range(n):
        result_y.append(b0 + b1 * plan_matrix[i][0] + b2 * plan_matrix[i][1] + b3 * plan_matrix[i][2] +
                        b12 * plan_matrix[i][3] + b13 * plan_matrix[i][4] + b23 * plan_matrix[i][5] +
                        b123 * plan_matrix[i][6])



def kohren():
    global m
    global gp
    global s
    s = [sum([(y_matrix[j][i] - average_y[i]) ** 2 for i in range(m)]) / m for j in range(n)]
    gp = max(s) / sum(s)
    khr = f.ppf(q=1 - q1, dfn=f1, dfd=(f2 - 1) * f1)
    gt = khr / (khr + f2 - 1)
    if gp > gt:
        m += 1
        main()


def studens():
    global n, m
    global average_y
    global b0, b1, b2, b3, b12, b13, b23, b123
    global d
    global blist
    global sb

    d = 8

    sb = sum(s) / n
    s_beta_2 = sb / (n * m)
    s_beta = s_beta_2 ** (1 / 2)

    bb = [b0, b1, b2, b3, b12, b13, b23, b123]
    t_tmp = [abs(bb[i]) / s_beta for i in range(n)]
    tt = t.ppf(q=0.975, df=f3)
    blist = [b0, b1, b2, b3, b12, b13, b23, b123]

    for i in range(n):
        if t_tmp[i] < tt:
            blist[i] = 0
            d -= 1

    return t_tmp


def fish():
    global n
    global d
    global blist
    global plan_matrix_for_output
    global average_y
    global sb
    global fp
    f4 = n - d
    y_reg = [b0 + b1 * plan_matrix[i][0] + b2 * plan_matrix[i][1] + b3 * plan_matrix[i][2] +
             b12 * plan_matrix[i][3] + b13 * plan_matrix[i][4] + b23 * plan_matrix[i][5] +
             b123 * plan_matrix[i][6] for i in range(n)]
    sad = (m / (n - d)) * int(sum([(y_reg[i] - average_y[i]) ** 2 for i in range(n)]))
    fp = sad / sb
    if fp > f.ppf(q=0.95, dfn=f4, dfd=f3):
        return 'The regression equation is inadequate to the original at a significance level of 0.05'
    else:
        return 'The regression equation is adequate to the original at a significance level of 0.05'


main()
kohren()
t_list = studens()
fisher = fish()




plan_table = BeautifulTable()
headers_x = ['X{}'.format(i) for i in range(0, m+1)]
headers_x.extend(['X12', 'X13', 'X23', 'X123'])
headers_y = ['Y{}'.format(i) for i in range(1, m+1)]
headers_y.extend(['av_Y', 'S^2'])
plan_table.columns.header = [*headers_x, *headers_y]
x0 = [[1] for _ in range(n)]
for i in range(n):
    plan_table.rows.append([*x0[i], *plan_matrix[i], *y_matrix[i], average_y[i], s[i]])
print('Матриця планування:')
print(plan_table)

norm_table = BeautifulTable()
headers_x = ['X{}'.format(i) for i in range(0, m+1)]
headers_x.extend(['X12', 'X13', 'X23', 'X123'])
headers_y = ['Y{}'.format(i) for i in range(1, m+1)]
headers_y.extend(['av_Y', 'S^2'])
norm_table.columns.header = [*headers_x, *headers_y]
x0 = [[1] for _ in range(n)]
for i in range(n):
    norm_table.rows.append([*x0[i], *plan_matrix_normal[i], *y_matrix[i], average_y[i], s[i]])

print('Нормована матриця:')
print(norm_table)

print('Kohren check')
print('Gp = {} < 0.7679'.format(round(gp, 3)))

print('Studens')
for i in range(len(t_list)):
    print('t{} = {}'.format(i, round(t_list[i], 3)))

print('Fisher')
print(fisher)

print('The equation')
print('y = {} + {} * x1 + {} * x2 + {} * x3 + {} * x1x2 + {} * x1x3 + {} * x2x3 + {} * x1x2x3'
      .format(round(b0, 3), round(b1, 3), round(b2, 3), round(b3, 3), round(b12, 3), round(b13, 3), round(b23, 3),
              round(b123, 3)))
