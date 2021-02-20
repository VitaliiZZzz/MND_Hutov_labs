import random

# Генерація масиву Х-ів
array_x = [[random.randint(0, 20) for i in range(3)] for j in range(8)]

# Генерація масиву змінних а
array_a = [random.randint(0, 20) for i in range(4)]

# Розрахунок У-ків
array_y = [array_a[0]+array_a[1]*array_x[i][0]+array_a[2]*array_x[i][1]+array_a[3]*array_x[i][2]
             for i in range(len(array_x))]

# Знаходження х0
array_average_x = [(min([(array_x[j][i]) for j in range(8)])+max([(array_x[j][i]) for j in range(8)]))/2 for i in range(3)]

# Обчислення інтервалу зміни dx
array_dx = [(array_average_x[i] - min([(array_x[j][i]) for j in range(8)])) for i in range(3)]

# Нормування Xn
array_xn = [[(array_x[j][i]-array_average_x[i])/array_dx[i] for i in range(3)] for j in range(8)]

# Знаходження точки плану
average_y = sum(array_y)/len(array_y)

difference_average_y = [array_y[i]-average_y for i in range(len(array_y))]
tmp_list = [min(difference_average_y) for i in range(len(array_y))]
for i in range(len(difference_average_y)):
    if difference_average_y[i] < 0:
        tmp_list[i] = difference_average_y[i]
index = None
for i in range(len(tmp_list)):
    if tmp_list[i] == max(tmp_list):
        index = i

print('_'*100)
print('Масив Х:')
for i in range(len(array_x)): print('[{:2}, {:2}, {:2}]'.format(array_x[i][0], array_x[i][1], array_x[i][2]))
print('_'*100)
print('a0: {}\na1: {}\na2: {}\na3: {}'.format(array_a[0], array_a[1], array_a[2], array_a[3]))
print('_'*100)
print('Масив Y:')
print(array_y)
print('_'*100)
print('Список значень x0:')
print(array_average_x)
print('_'*100)
print('Список значень dx:')
print(array_dx)
print('_'*100)
print('Нормовані значення Хn:')
for i in range(len(array_xn)): print('[{:7.4f}, {:7.4f}, {:7.4f}]'.format(array_xn[i][0], array_xn[i][1], array_xn[i][2]))
print('_'*100)
print('Середнє значення Y:', average_y)
print('_'*100)
print('Для →Y - де Y середнє значення масиву Y')
print('x0: {}\nx1: {}\nx2: {}'.format(array_x[index][0], array_x[index][1], array_x[index][2]))
print('Y = ', array_y[index])
print('_'*100)
