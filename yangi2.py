import numpy as np
import math
import matplotlib.pyplot as plt

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


def find_f(y_l):
    return (1 / 8) * (y_l + 1)


find_f_massiv = np.zeros(N + 1)

for i in range(0, N + 1):
    find_f_massiv[i] = find_f(y_l_array[i])


def find_b_n(n):
    if (n == 0 or n == N):
        c_n = 2
    else:
        c_n = 1

    sum = 0
    for i in range(0, N + 1):
        if (i == 0 or i == N):
            c_l = 2
        else:
            c_l = 1
        sum += (1 / c_l) * find_f_massiv[i] * find_T_n_matrix[n, i]

    return 2 / (N * c_n) * sum


a_matrix = np.zeros((N + 1, N + 1))
free_b_array = np.zeros(N + 1)

for i in range(0, N - 1):
    if i == 0:
        c_i = 2
    else:
        c_i = 1

    for p in range(i + 2, N + 1):
        if (p - i) % 2 == 0:
            a_matrix[i, p] += 1 / c_i * p * (pow(p, 2) - pow(i, 2))
    a_matrix[i, i] -= 1 / (2 * c_i)
    free_b_array[i] = find_b_n(i)

for i in range(0, N + 1):
    a_matrix[N - 1, i] = 1

    if i % 2 == 0:
        a_matrix[N, i] = 1
    else:
        a_matrix[N, i] = -1

print("a_matrix")
print(a_matrix)
print("free b bektor")
print(free_b_array)

a_gauss_result = np.linalg.solve(a_matrix, free_b_array)
print("gauss yechim")
print(a_gauss_result)


def find_U_taqribiy(y_l):
    sum = 0
    for i in range(0, N + 1):
        sum += a_gauss_result[i] * find_T_n_matrix[i, y_l]
    return sum


find_U_taqribiy_array = np.zeros(N + 1)

for i in range(N + 1):
    find_U_taqribiy_array[i] = find_U_taqribiy(i)


def find_U_aniq_yechim(y_l):
    return ((eps - 0.5) / (1 - math.pow(math.e, ((-1 / eps)))) *
            (1 - math.pow(math.e, (-(y_l_array[y_l] + 1) / (2 * eps)))) - eps * ((y_l_array[y_l] + 1) / 2) +
            (math.pow((y_l_array[y_l] + 1), 2) / 8))


find_U_aniq_yechim_massiv = np.zeros(N + 1)
for i in range(N + 1):
    find_U_aniq_yechim_massiv[i] = find_U_aniq_yechim(i)

find_xatolik_massiv = np.zeros(N + 1)

for i in range(N + 1):
    find_xatolik_massiv[i] = math.fabs(find_U_aniq_yechim_massiv[i] - find_U_taqribiy_array[i])

print(find_xatolik_massiv)

plt.plot(find_U_aniq_yechim_massiv, 'o-', label="Aniq yechim", color="green")
plt.plot(find_U_taqribiy_array, 'o--', label="Taqribiy yechim", color="red")
plt.title("Aniq va Taqribiy yechimlar grafigi")
plt.xlabel(f"N = {N},  eps = {eps}", loc="right", color="red", size=15)
plt.ylabel("natijalar")
plt.legend(loc="best", fontsize=12)
plt.grid(color="blue", linestyle=":", linewidth=0.5)

plt.show()
