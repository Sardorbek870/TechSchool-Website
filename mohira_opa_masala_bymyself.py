import numpy as np
import math
pi = math.pi
e = np.e
N = int(input("N = "))
eps = float(input("eps = "))

l_values = np.arange(N + 1)
y_l_array = -1 * np.cos(pi * l_values / N)

print("y(l) larning qiymati")
for i in range(0, N + 1):
    print("y(", i, ")=", y_l_array[i])
print("-------------------------------------------------")


find_T_n_matrix = np.zeros((N + 1, N + 1))
find_T_n_matrix[0, :] = 1
find_T_n_matrix[1, :] = y_l_array
for n in range(2, N + 1):
    find_T_n_matrix[n, :] = 2 * y_l_array * find_T_n_matrix[n - 1, :] - find_T_n_matrix[n - 2, :]

print("T larning qiymatlari")
for i in range(0, N + 1):
    for j in range(0, N + 1):
        print("T", i, "(y(", j, ")) = ", find_T_n_matrix[i, j])
    print("")
print("---------------------------------------------------")


def find_f_y_l(y_l):
    return eps * ((pow(eps,4) + pow(eps, 2) + eps) * pow(y_l_array[y_l], 4)
                  + (16 * pow(eps,3) + 8 * eps + 4) * pow(y_l_array[y_l], 3)
                  + (70 * pow(eps, 2) - 2 * pow(eps, 4) - eps + 12) * pow(y_l_array[y_l], 2)
                  + (88 * eps - 16 * pow( eps, 3) - 4) * y_l_array[y_l]
                  + (pow(eps, 4) - 23 * pow(eps, 2) + 20)) * math.exp(eps * y_l_array[y_l])


find_f_y_l_array = np.zeros(N+1)
for i in range(N+1):
    find_f_y_l_array[i] = find_f_y_l(i)

for i in range(0, N + 1):
    print("f(y(", i, ")) = ", find_f_y_l_array[i])
print("----------------------------------------------------")


def find_b_n(n):
    if n == 0 or n == N:
        c_n = 2
    else:
        c_n = 1

    sum = 0
    for i in range(0, N + 1):
        if(i == 0 or i == N):
            c_l = 2
        else:
            c_l = 1
        sum += (2 / c_l) * find_f_y_l_array[i] * find_T_n_matrix[n, i]

    return 2 / (N * c_n) * sum


find_b_n_array = np.zeros(N + 1)
for i in range(N + 1):
    find_b_n_array[i] = find_b_n(i)


a_matrix = np.zeros((N + 1, N + 1))
free_b_array = np.zeros(N + 1)

for n in range(0, N - 3):
    if n == 0:
        c_n = 2
    else:
        c_n = 1

    for p in range(n+4, N + 1):
        if p % 2 == n:
            a_matrix[n][p] += (eps / (24 * c_n) * p *
                               (p**2 * (pow((p ** 2 - 4), 2) - 3 * n ** 2 * pow(p,4) +
                                        3 * pow(n,4) * p**2 - n ** 2 * pow((n**2 - 4), 2))))

    for p in range(n+2, N+1):
        if p % 2 == n:
            a_matrix[n][p] -= (2/c_n * p * (p**2 - n**2))

    a_matrix[n][n] += 1
free_b_array[n] += find_b_n_array[n]

for n in range(0, N+1):
    if n % 2 == 0:
        a_matrix[N - 3][n] = 1
        a_matrix[N - 2][n] = n * n
    if n % 2 == 1:
        a_matrix[N - 1][n] = 1
        a_matrix[N][n] = n * n

free_b_array[N - 3] = 0
free_b_array[N - 2] = 0
free_b_array[N - 1] = 0
free_b_array[N] = 0

print("a matritsaning koeffitsientlari")
print(a_matrix)
print("-----------------------------------------------")
print("b ozod hadlar vektori")
print(free_b_array)
print("-------------------------------------------")


def gauss_method(a_matrix, free_b_array):
    return np.linalg.solve(a_matrix, free_b_array)

x_gauss_yechim = gauss_method(a_matrix, free_b_array)

print("gauss usuli orqali olingan yechimlar ")
print(x_gauss_yechim)
print("-----------------------------------------------")

a_n_1 = np.zeros(N + 1)
for n in range(0, N+1):
    if n == 0 or n == N:
        c_n = 2
    else:
        c_n = 1
    for p in range(n+1, N+1):
        if (p + n - 1) % 2 == 0:
            a_n_1 += 2/ c_n * p * x_gauss_yechim[p]


hosila_1 = np.zeros(N+1)
for i in range(N+1):
    natija = 0
    for j in range(0, N+1):
        natija += a_n_1[j] * find_T_n_matrix[j,i]
    hosila_1[i] = natija

print("hosila 1")
print(hosila_1)
print("------------------------------------------------------")

a_n_2 = np.zeros(N + 1)

for n in range(0, N+1):
    if n == 0 or n == N:
        c_n = 2
    else:
        c_n = 1
    for p in range(n+2, N+1):

        if (p - n)  % 2 == 0:
            a_n_1 += 1 / c_n * p * (p**2 - n**2) * x_gauss_yechim[p]
hosila_2 = np.zeros(N+1)

for i in range(N+1):
    natija = 0
    for j in range(N+1):
        natija += a_n_2[j] * find_T_n_matrix[j,i]
    hosila_2[i] = natija

print("hosila 2")
print(hosila_2)
print("------------------------------------------------------")

a_n_3 = np.zeros(N + 1)
for n in range(0, N+1):
    if n == 0 or n == N:
        c_n = 2
    else:
        c_n = 1
    for p in range(n+3, N+1):
        if (p + n) % 2 == 1:
            a_n_3 += 1/ (4 * c_n) * p * ((p**2 - 1)**2  - 2 * (p**2 +1) * n**2 + pow(n,4)) * x_gauss_yechim[p]

hosila_3 = np.zeros(N+1)

for i in range(N+1):
    natija = 0
    for j in range(N+1):
        natija += a_n_3[j] * find_T_n_matrix[j,i]
    hosila_3[i] = natija

print("hosila 3")
print(hosila_3)
print("------------------------------------------------------")

a_n_4 = np.zeros(N + 1)
for n in range(0, N+1):
    if n == 0 or n == N:
        c_n = 2
    else:
        c_n = 1
    for p in range(n+4, N+1):
        if p % 2 == n:
            a_n_1 += 1/(24* c_n) * p * (pow(p,2) * pow((p**2 - 4),2) -
                                        3 * n**2 * pow(p,4) + 3 * n**4 * p**2 - n**2 *
                                        pow((n**2 - 4),2)) * x_gauss_yechim[p]

hosila_4 = np.zeros(N+1)

for i in range(N+1):
    natija = 0
    for j in range(N+1):
        natija += a_n_4[j] * find_T_n_matrix[j,i]
    hosila_4[i] = natija


print("hosila 4")
print(hosila_4)
print("------------------------------------------------------")
aniq_yechim = np.zeros(N+1)
hosila_aniq_yechim_1 = np.zeros(N+1)
hosila_aniq_yechim_2 = np.zeros(N+1)
hosila_aniq_yechim_3 = np.zeros(N+1)
hosila_aniq_yechim_4 = np.zeros(N+1)

for i in range(N+1):
    aniq_yechim[i] = (1 - pow(y_l_array[i], 2)) * pow(e, eps * y_l_array[i])
    hosila_aniq_yechim_1[i] = (-eps * pow(y_l_array[i], 2) - 2 * y_l_array[i] + eps) * pow(e, eps * y_l_array[i])
    hosila_aniq_yechim_2[i] = (-1 * pow(eps, 2) * pow(y_l_array[i], 2) - 4 * eps * y_l_array[i] + pow(eps, 2) - 2) * pow(e, eps * y_l_array[i])
    hosila_aniq_yechim_3[i] = -eps * (pow(eps, 2) * pow(y_l_array[i], 2) + 6 * eps * y_l_array[i] - pow(eps, 2) + 6) * pow(e, eps * y_l_array[i])
    hosila_aniq_yechim_4[i] = -1 * pow(eps, 2) * (pow(eps, 2) * pow(y_l_array[i], 2) + 8 * eps * y_l_array[i] - pow(eps, 2) + 12) * pow(e, eps * y_l_array[i])

print("1 - tartibli xosila dagi hotoliklar ")
for i in range(N+1):
    print(f"hosila1[{i}] = {hosila_1[i]}    -    aniq_hosila1[{i}] = {hosila_aniq_yechim_1[i]},   xatolik = {np.fabs(hosila_aniq_yechim_1[i] - hosila_1[i])}")
print("------------------------------------")

print("2 - tartibli xosila dagi hotoliklar ")
for i in range(N+1):
    print(f"hosila2[{i}] = {hosila_2[i]}    -    aniq_hosila2[{i}] = {hosila_aniq_yechim_2[i]},   xatolik = {np.fabs(hosila_aniq_yechim_2[i] - hosila_2[i])}")
print("------------------------------------")

print("3 - tartibli xosila dagi hotoliklar ")
for i in range(N+1):
    print(f"hosila3[{i}] = {hosila_3[i]}    -    aniq_hosila3[{i}] = {hosila_aniq_yechim_3[i]},   xatolik = {np.fabs(hosila_aniq_yechim_3[i] - hosila_3[i])}")
print("------------------------------------")

print("4 - tartibli xosila dagi hotoliklar ")
for i in range(N+1):
    print(f"hosila4[{i}] = {hosila_4[i]}    -    aniq_hosila4[{i}] = {hosila_aniq_yechim_4[i]},   xatolik = {np.fabs(hosila_aniq_yechim_4[i] - hosila_4[i])}")
print("------------------------------------")



