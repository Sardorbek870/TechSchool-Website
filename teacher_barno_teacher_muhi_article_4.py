import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
e = np.e

N = int(input("N = "))
eps = float(input("eps = "))

l_values = np.arange(N+1)
y_l_array = -1 * np.cos(pi * l_values / N)

T_n_matrix = np.zeros((N+1, N+1))
T_n_matrix[0,:] = 1
T_n_matrix[1, :] = y_l_array
for n in range(2, N+1):
    T_n_matrix[n, :] = 2 * y_l_array * T_n_matrix[n-1, :] - T_n_matrix[n-2, :]

U_y_aniq_yechimlar = np.zeros(N+1)
for i in range(N+1):
    U_y_aniq_yechimlar[i] = pow((1 - pow(y_l_array[i], 2)),2) * np.exp(eps * y_l_array[i])

for i in range(N+1):
    print(f"{i}) {y_l_array[i]} -> {U_y_aniq_yechimlar[i]}")

import math
def f(y):
    return ((eps ** 5 - 2 * eps ** 2) * y ** 4 + (16 * eps ** 4 - 16 * eps) * y ** 3 + (
                72 * eps ** 3 - 2 * eps ** 5 - 24 + 4 * eps ** 2) * y ** 2 + (
                        96 * eps ** 2 - 16 * eps ** 4 + 16 * eps) * y - 2 * eps ** 2 + 8) * math.exp(eps * y)


f_y_array = np.zeros(N+1)
for i in range(N+1):
    #f_y_array[i] = eps * ((pow(eps,4) + pow(eps,2) + eps) * pow(y_l_array[i],4) + (16 * pow(eps,3) + 8 * eps + 4) * pow(y_l_array[i],3) + (70 * pow(eps,2) - 2*pow(eps,4) - eps + 12) * pow(y_l_array[i],2) + (88 * eps - 16 * pow(eps, 3) - 4) * y_l_array[i] + pow(eps,4) - 23 * pow(eps,2) + 20) * np.exp(eps * y_l_array[i])
    f_y_array[i] = f(y_l_array[i])



c_array = np.ones(N+1)
c_array[0] = 2
c_array[N] = 2

def find_b_n(n):
    sum = 0
    for i in range(N+1):
        sum += (1/c_array[i]) * f_y_array[i] * T_n_matrix[n,i]

    return 2 / (N * c_array[n]) * sum


a_matrix = np.zeros((N+1, N+1))
free_b_array = np.zeros(N+1)

for n in range(N-3):
    for p in range(n+4, N+1):
        if (p-n) % 2 == 0:
            a_matrix[n,p] += eps / (24 * c_array[n]) * p * (pow(p,2) * pow((p*p - 4),2) - 3*n*n*pow(p,4) + 3*pow(n,4)*p*p - n*n*pow((n*n - 4),2))

    for p in range(n+2, N+1):
        if (p-n) % 2 == 0:
            a_matrix[n,p] -= (2 / c_array[n]) * p * (p*p - n*n)

    a_matrix[n,n] += 1
    free_b_array[n] = find_b_n(n)

for n in range(N+1):
    if n % 2 == 0:
        a_matrix[N-3,n] = 1
        a_matrix[N-2, n] = pow(n,2)

for n in range(1, N+1):

    if (n - 1) % 2 == 0:
        a_matrix[N - 1, n] = 1
        a_matrix[N, n] = pow(n,2)


a_gauss_result = np.linalg.solve(a_matrix, free_b_array)

def find_U_taqribiy(y_l):
    sum = 0
    for i in range(0, N+1):
        sum += a_gauss_result[i] * T_n_matrix[i, y_l]
    return sum

find_U_taqribiy_array = np.zeros(N+1)

for i in range(N+1):
    find_U_taqribiy_array[i] = find_U_taqribiy(i)

print("Taqribiy yechimlar")
for i in range(N+1):
    print(f"U(y({i})) = {find_U_taqribiy_array[i]}")
print("-" * 20)

plt.plot(y_l_array, U_y_aniq_yechimlar, 'o--', label="Aniq yechim", color="green")
plt.plot(y_l_array, find_U_taqribiy_array, 'o--', label="Taqribiy yechim", color="red")
plt.title("Aniq va Taqribiy yechimlar grafigi")
plt.xlabel(f"N = {N},  eps = {eps}", loc="right", color="red", size=15)
plt.ylabel("natijalar")
plt.legend(loc="best", fontsize=12)
plt.grid(color="blue", linestyle=":", linewidth=0.5)

plt.show()