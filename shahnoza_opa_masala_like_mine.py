# c_N = 1 va b0 = 1/2 * b0 ga o'zgartirilgan holat
import math
import sys
import numpy
import matplotlib.pyplot as plt
import numpy as np


N = int(input("N = "))

eps = float(input("eps = "))
pi = math.pi


def find_y_l(l):
    return -1 * math.cos(pi * l / N)
    #return math.cos(pi *l / N)

find_y_l_matrix = numpy.zeros(N+1)
for i in range(0,N+1):
    find_y_l_matrix[i] = find_y_l(i)


def find_T_n(n, y_l):
    if n == 0:
        return 1
    if (n == 1):
        return y_l
    result = 2 * y_l * find_T_n(n - 1, y_l) - find_T_n(n - 2, y_l)
    return result

find_T_n_massiv = numpy.zeros((N+1, N+1))
for i in range(0, N+1):
    for j in range(N+1):
        find_T_n_massiv[i, j] = find_T_n(i, find_y_l_matrix[j])


l_values = np.arange(N + 1)
find_y_l_matrix = -1 * np.cos(pi * l_values / N)


n_values = np.arange(N + 1)
find_T_n_massiv = np.zeros((N + 1, N + 1))

find_T_n_massiv[0, :] = 1
find_T_n_massiv[1, :] = find_y_l_matrix

for n in range(2, N + 1):
    find_T_n_massiv[n, :] = 2 * find_y_l_matrix * find_T_n_massiv[n - 1, :] - find_T_n_massiv[n - 2, :]


def find_f(y_l):
    return (1 / 8) * (y_l + 1)

find_f_massiv = numpy.zeros(N+1)

for i in range(0,N+1):
    find_f_massiv[i] = find_f(find_y_l_matrix[i])



def find_b_n(n):
    if (n == 0):
        c_n = 2
    else:
        c_n = 1

    sum = 0
    for i in range(0, N + 1):
        if (i == 0 or i == N) :
            c_l = 2
        else:
            c_l = 1
        sum += (1 / c_l) * find_f_massiv[i] * find_T_n_massiv[n, i]

    return 2 / (N * c_n) * sum


find_b_n_matrix = ([0 for i in range(0,N+1)])
for i in range(0,N+1):
    find_b_n_matrix[i] = find_b_n(i)



a_matrix = numpy.zeros((N + 1, N + 1))
free_b_vector = numpy.zeros(N + 1)

for i in range(0, N + 1):
    if (i == 0):
        a_matrix[0][i] = 1 / 2   # mistake here
       # a_matrix[0][i] = 1
    elif (i % 2 == 0):
        a_matrix[0][i] = 1
    else:
        a_matrix[1][i] = 1

free_b_vector[0] = 0
free_b_vector[1] = 0

for i in range(2, N + 1):

    if i == N - 1 or i == N:
        free_b_vector[i] =  (i + 1) * find_b_n_matrix[i - 2] - 2 * i * find_b_n_matrix[i]
    else:
        free_b_vector[i] = ((i + 1) * find_b_n_matrix[i - 2] -
                             2 * i * find_b_n_matrix[i] + (i - 1) *
                            find_b_n_matrix[i + 2])



    a_matrix[i][i - 1] = i * i - 1
    a_matrix[i][i] = 4 * i * (i * i - 1) * eps

    lastIndex = i + 1
    if (lastIndex <= N):
        a_matrix[i][lastIndex] = -1 * (i * i - 1)


def gauss_method(a_matrix, free_b_matrix):
    return numpy.linalg.solve(a_matrix, free_b_matrix)
    # column_vector = free_b_matrix.reshape(-1, 1)
    # # Add the column vector to the array
    # a = numpy.hstack((a_matrix, column_vector))
    #
    # x = numpy.zeros(N + 1)
    # n = N + 1
    # # Applying Gauss Elimination
    # for i in range(n):
    #     if a[i][i] == 0.0:
    #         sys.exit('Divide by zero detected!')
    #
    #     for j in range(i + 1, n):
    #         ratio = a[j][i] / a[i][i]
    #
    #         for k in range(n + 1):
    #             a[j][k] = a[j][k] - ratio * a[i][k]
    #
    # # Back Substitution
    # x[n - 1] = a[n - 1][n] / a[n - 1][n - 1]
    #
    # for i in range(n - 2, -1, -1):
    #     x[i] = a[i][n]
    #
    #     for j in range(i + 1, n):
    #         x[i] = x[i] - a[i][j] * x[j]
    #
    #     x[i] = x[i] / a[i][i]
    #
    # return x


x_gauss_yechim = gauss_method(a_matrix, free_b_vector)


def find_U_taqribiy(y_l):
    sum = 0
    for i in range(0, N + 1):
        if (i == 0):
            sum += 1 / 2 * x_gauss_yechim[i] * find_T_n_massiv[i, y_l]
        else:
            sum += x_gauss_yechim[i] * find_T_n_massiv[i, y_l]

    return sum

find_U_taqribiy_massiv = numpy.zeros(N+1)

for i in range(N+1):
    find_U_taqribiy_massiv[i] = find_U_taqribiy(i)


def find_U_aniq_yechim(y_l):
    return ((eps - 0.5) / (1 - math.pow(math.e, ((-1 / eps)))) *
            (1 - math.pow(math.e, (-(find_y_l_matrix[y_l] + 1) / (2 * eps)))) - eps * ((find_y_l_matrix[y_l] + 1) / 2) +
            (math.pow((find_y_l_matrix[y_l] + 1), 2) / 8))


find_U_aniq_yechim_massiv = numpy.zeros(N+1)
for i in range(N+1):
    find_U_aniq_yechim_massiv[i] = find_U_aniq_yechim(i)





# def find_xatolik(aniq_yechim, taqribiy_yechim):
#     return math.fabs(aniq_yechim - taqribiy_yechim)


find_xatolik_massiv = numpy.zeros(N+1)

for i in range(N+1):
    find_xatolik_massiv[i] = math.fabs(find_U_aniq_yechim_massiv[i] - find_U_taqribiy_massiv[i])



print("y(l) larning qiymati")
for i in range(0, N + 1):
    print("y(", i, ")=", find_y_l_matrix[i])
print("-------------------------------------------------")
print("T larning qiymatlari")
for i in range(0, N + 1):
    for j in range(0, N + 1):
        print("T", i, "(y(", j, ")) = ", find_T_n_massiv[i, j])

    print("")
print("---------------------------------------------------")

for i in range(0, N + 1):
    print("f(y(", i, ")) = ", find_f_massiv[i])
print("----------------------------------------------------")


for i in range(0, N + 1):
    print("b(", i, ") = ", find_b_n_matrix[i])

print("--------------------------------------------")
print("a matritsaning koeffitsientlari")
print(a_matrix)
print("-----------------------------------------------")
print("b ozod hadlar vektori")
print(free_b_vector)
print("-------------------------------------------")

solution = gauss_method(a_matrix, free_b_vector)
print("gauss usuli orqali olingan yechimlar ")
print(solution)
print("-----------------------------------------------")
print("taqribiy yechim")
for i in range(0, N + 1):
    print("u(y(", i, ")) = ", find_U_taqribiy_massiv[i])
print("-----------------------------------------------")

print("aniq yechimlar")
for i in range(0, N + 1):
    print("u(y(", i, ")) = ", find_U_aniq_yechim_massiv[i])
print("--------------------------------------------------")



print("xatolik")
for i in range(0, N + 1):

    print("[", find_U_aniq_yechim_massiv[i], " - ", find_U_taqribiy_massiv[i], "] = ", find_xatolik_massiv[i])

print("##############################################################################!")
for i in range(N):
    natija = 0
    for j in range(N):
        natija += find_b_n_matrix[j] * find_T_n_massiv[j, i]
    print(f"f(y{i}) = {natija}  -  f(y{i}) = {find_f_massiv[i]} =  {natija - find_f_massiv[i]}")
print("###############################################################################")

plt.plot(find_U_aniq_yechim_massiv, 'o-',  label = "Aniq yechim", color = "green")
plt.plot(find_U_taqribiy_massiv,'o--', label = "Taqribiy yechim", color = "red")
plt.title("Aniq va Taqribiy yechimlar grafigi")
plt.xlabel(f"N = {N},  eps = {eps}", loc = "right", color = "red", size = 15)
plt.ylabel("natijalar")
plt.legend(loc = "best", fontsize = 12)
plt.grid(color = "blue", linestyle = ":", linewidth = 0.5)

plt.show()