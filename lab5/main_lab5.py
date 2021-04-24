import random
import sklearn.linear_model as lm
from scipy.stats import f, t
from functools import partial
from pyDOE2 import *
from beautifultable import BeautifulTable

# Гутов Віталій
# Варіант 108:
# x1_min = -5, x1_max = 7,
# x2_min = -10, x2_max = 3,
# x3_min = -7, x3_max = 1
# y_min = 200 + xc_min
# y_max = 200 + xc_max


def plan_matrix(n, m):
    y = np.zeros(shape=(n, m))
    for i in range(n):
        for j in range(m):
            y[i][j] = random.randint(y_min, y_max)

    if n > 14:
        no = n - 14
    else:
        no = 1
    x_norm = ccdesign(3, center=(0, no))
    x_norm = np.insert(x_norm, 0, 1, axis=1)

    for i in range(4, 11):
        x_norm = np.insert(x_norm, i, 0, axis=1)

    l = 1.215
    for i in range(len(x_norm)):
        for j in range(len(x_norm[i])):
            if x_norm[i][j] < -1 or x_norm[i][j] > 1:
                if x_norm[i][j] < 0:
                    x_norm[i][j] = -l
                else:
                    x_norm[i][j] = l

    def add_sq_nums(x):
        for i in range(len(x)):
            x[i][4] = x[i][1] * x[i][2]
            x[i][5] = x[i][1] * x[i][3]
            x[i][6] = x[i][2] * x[i][3]
            x[i][7] = x[i][1] * x[i][3] * x[i][2]
            x[i][8] = x[i][1] ** 2
            x[i][9] = x[i][2] ** 2
            x[i][10] = x[i][3] ** 2
        return x

    x_norm = add_sq_nums(x_norm)
    x = np.ones(shape=(len(x_norm), len(x_norm[0])), dtype=np.int64)
    for i in range(8):
        for j in range(1, 4):
            if x_norm[i][j] == -1:
                x[i][j] = x_range[j - 1][0]
            else:
                x[i][j] = x_range[j - 1][1]

    for i in range(8, len(x)):
        for j in range(1, 3):
            x[i][j] = (x_range[j - 1][0] + x_range[j - 1][1]) / 2
    dx = [x_range[i][1] - (x_range[i][0] + x_range[i][1]) / 2 for i in range(3)]
    x[8][1] = l * dx[0] + x[9][1]
    x[9][1] = -l * dx[0] + x[9][1]
    x[10][2] = l * dx[1] + x[9][2]
    x[11][2] = -l * dx[1] + x[9][2]
    x[12][3] = l * dx[2] + x[9][3]
    x[13][3] = -l * dx[2] + x[9][3]
    x = add_sq_nums(x)
    x_table = BeautifulTable()
    for i in range(n):
        x_table.rows.append([*x[i]])
    print('x matrix:')
    print(x_table)
    x_norm_table = BeautifulTable()
    for i in range(n):
        x_norm_table.rows.append([*x_norm[i]])
    print('Normalized x matrix:')
    print(x_norm_table)
    return x, y, x_norm


def regression(x, b):
    y = sum([x[i] * b[i] for i in range(len(x))])
    return y


def s_kv(y, y_aver, n, m):
    res = []
    for i in range(n):
        s = sum([(y_aver[i] - y[i][j]) ** 2 for j in range(m)]) / m
        res.append(round(s, 3))
    return res


def coef_finding(x, y, norm=False):
    skm = lm.LinearRegression(fit_intercept=False)
    skm.fit(x, y)
    b = skm.coef_
    if norm == 1:
        print('\nCoefficients of the regression equation with normalized x:')
    else:
        print('\nCoefficients of the regression equation:')
    b = [round(i, 3) for i in b]
    print(b)
    print('\nThe result of the equation with the found coefficients:\n{}'.format(np.dot(x, b)))
    return b


def kohren_kr(y, y_aver, n, m):
    f1 = m - 1
    f2 = n
    q = 0.05
    skv = s_kv(y, y_aver, n, m)
    gp = max(skv) / sum(skv)
    print('\nKohren check')
    return gp


def kohren(f1, f2, q=0.05):
    q1 = q / f1
    fisher_value = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)
    return fisher_value / (fisher_value + f1 - 1)


def bs(x, y_aver, n):
    res = [sum(1 * y for y in y_aver) / n]
    for i in range(len(x[0])):
        b = sum(j[0] * j[1] for j in zip(x[:, i], y_aver)) / n
        res.append(b)
    return res


def student_kr(x, y, y_aver, n, m):
    skv = s_kv(y, y_aver, n, m)
    skv_aver = sum(skv) / n
    sbs_tmp = (skv_aver / n / m) ** 0.5
    bs_tmp = bs(x, y_aver, n)
    ts = [round(abs(b) / sbs_tmp, 3) for b in bs_tmp]
    return ts


def fisher_kr(y, y_aver, y_new, n, m, d):
    S_ad = m / (n - d) * sum([(y_new[i] - y_aver[i]) ** 2 for i in range(len(y))])
    skv = s_kv(y, y_aver, n, m)
    skv_aver = sum(skv) / n
    return S_ad / skv_aver


def check(x, y, b, n, m):
    print('\nCheck the equation:')
    f1 = m - 1
    f2 = n
    f3 = f1 * f2
    q = 0.05
    student = partial(t.ppf, q=1 - q)
    t_student = student(df=f3)
    g_kr = kohren(f1, f2)
    y_aver = [round(sum(i) / len(i), 3) for i in y]
    print('\nThe average value of y:', y_aver)
    disp = s_kv(y, y_aver, n, m)
    print('The variance y:', disp)
    gp = kohren_kr(y, y_aver, n, m)
    print(f'gp = {gp}')
    if gp < g_kr:
        print('With a probability of {} dispersions are homogeneous.'.format(1 - q))
    else:
        print("It is necessary to increase the number of experiments")
        m += 1
        main(n, m)
    ts = student_kr(x[:, 1:], y, y_aver, n, m)
    print('\nStudent criterion:\n{}:'.format(ts))
    res = [t for t in ts if t > t_student]
    final_k = [b[i] for i in range(len(ts)) if ts[i] in res]
    print('\nThe coefficients {} are statistically insignificant, so we exclude them from the equation.'.format([round(i, 3) for i in b if i not in final_k]))
    y_new = []
    for j in range(n):
        y_new.append(round(regression([x[j][i] for i in range(len(ts)) if ts[i] in res], final_k), 3))
    print('The value of y with coefficients {}: '.format(final_k))
    print(y_new)
    d = len(res)
    if d >= n:
        print('\nF4 <= 0')
        print('')
        return
    f4 = n - d
    f_p = fisher_kr(y, y_aver, y_new, n, m, d)
    fisher = partial(f.ppf, q=0.95)
    f_t = fisher(dfn=f4, dfd=f3)
    print('\nFisher adequacy check')
    print('fp =', f_p)
    print('ft =', f_t)
    if f_p < f_t:
        print('The mathematical model is adequate to the experimental data')
    else:
        print('The mathematical model is inadequate to the experimental data')


def main(n, m):
    x, y, x_norm = plan_matrix(n, m)
    y5_aver = [round(sum(i) / len(i), 3) for i in y]
    b = coef_finding(x, y5_aver)
    check(x_norm, y, b, n, m)


x_range = ((-5, 7), (-10, 3), (-7, 1))
x_aver_max = sum([x[1] for x in x_range]) / 3
x_aver_min = sum([x[0] for x in x_range]) / 3
y_min = 200 + int(x_aver_min)
y_max = 200 + int(x_aver_max)


main(15, 3)
