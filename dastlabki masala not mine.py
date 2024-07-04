
import math
import numpy as np
import matplotlib.pyplot as plt
N = int(input("N = "))
A = float(input("A = "))
#A = 0.5
eps = 0.01
pi = math.pi
e = math.e
def find_X_l(l):
    return -math.cos(pi * l / N)
find_X_l_array = np.zeros(N+1)
for i in range(N+1):
    find_X_l_array[i] = find_X_l(i)

def find_aniqYechim_U_x(x):
    return  (1 - x * x) * math.pow(math.e, (A * x))

find_aniqyechim_array = np.zeros(N + 1)
for i in range(N+1):
    find_aniqyechim_array[i] = find_aniqYechim_U_x(find_X_l_array[i])


# def find_T_X_l(n, x_l):#     if n == 0:
#         return 0#     if n == 1:
#         return x_l#     #return 2 * x_l * find_T_X_l(n - 1, x_l) - find_T_X_l(n - 2, x_l)
find_T_X_l_array = np.zeros((N + 1, N + 1))

find_T_X_l_array[0, :] = 1
find_T_X_l_array[1, :] = find_X_l_array
for n in range(2, N+1):
    find_T_X_l_array[n, :] = 2 * find_X_l_array * find_T_X_l_array[n - 1, :] - find_T_X_l_array[n - 2, :]
def find_f_x_l(x_l):
    # return (-1 *    #         ((-2 * (x_l * A + 1) + ((1 - x_l * x_l)
    #         * A - 2 * x_l)) * A)    #         * math.exp(A * x_l))
    #return -1 * (-2 * (x_l * A + 1) + ((1 - x_l * x_l) * A - 2 * x_l) * A) * math.pow(e, A * x_l)
    return (2 * (1 + 2 * A * x_l) - (1 - x_l * x_l) * A * A) * pow(2,(A * x_l))
    #return (pow(A,2) * pow(x_l, 2) + (4*A - 2*pow(A,2)) * x_l + pow(A,2) - 4 * A + 2) * np.exp(A * x_l)

find_f_x_l_array = np.zeros(N+1)
for i in range(N+1):
    find_f_x_l_array[i] = find_f_x_l(find_X_l_array[i])

def find_b_n(n):
    if n == 0 or n == N:
        c_k = 2
    else:
        c_k = 1
    sum = 0
    for i in range(0, N+1):
        if i == 0 or i == N:
            c_l = 2
        else:
            c_l = 1
        sum += 1 / c_l * find_f_x_l_array[i] * find_T_X_l_array[n, i]
    return 2 / (N * c_k) * sum


find_b_n_array = np.zeros(N+1)
for i in range(N+1):
    find_b_n_array[i] = find_b_n(i)


a_matrix = np.zeros((N+1, N+1))
for i in range(N+1):
    if i == 0:
        a_matrix[0][i] = 1 / 2
    elif i % 2 == 0:
        a_matrix[0][i] = 1
    else:
        a_matrix[1][i] = 1

b_vektor = np.zeros(N+1)
b_vektor[0] = 0
b_vektor[1] = 0
for i in range(2, N+1):
    a_matrix[i][i] = 1
    b_vektor[i] = (-1 / (4 * i * (i * i - 1))) * ((i + 1) *    find_b_n_array[i - 2] - 2 * i * find_b_n_array[i] +
    (i - 1) * (find_b_n_array[i+2] if not(i == N or i == N - 1) else 0))
def gauss_method(a_matrix, b_vektor):
    return np.linalg.solve(a_matrix, b_vektor)
a_gauss_yechim_vektor = gauss_method(a_matrix,b_vektor)
def find_U_x_l_taqribiy_yechim(x_l):
    sum = 0
    for i in range(N+1):
        sum += (1 / 2 if i == 0 else 1) * a_gauss_yechim_vektor[i] * find_T_X_l_array[i, x_l]
    return sum


find_U_x_l_taqribiy_yechim_array = np.zeros(N+1)

for i in range(N+1):
    find_U_x_l_taqribiy_yechim_array[i] = find_U_x_l_taqribiy_yechim(i)
find_xatoliklar_array = np.zeros(N+1)
for i in range(N+1):
    find_xatoliklar_array[i] = math.fabs(find_U_x_l_taqribiy_yechim_array[i] - find_aniqyechim_array[i])


print("x(l) larnig qiymatlari")
print(find_X_l_array)
print("---------------------------------------------")
print("Tn(x(l)) ning qiymatlari")
print(find_T_X_l_array)
print(("-------------------------------------------"))
print("b(k) ning qiymatlari")
print(find_b_n_array)
print("----------------------------------------------")
print("a matritsa ")
print(a_matrix)
print("---------------------------------------------")
print("b ozod hadlar vektori")
print(b_vektor)
print("---------------------------------------------")
print("gauss dan olingan yechim")
print(a_gauss_yechim_vektor)
print("---------------------------------------------")
print("aniq yechimlar")
print(find_aniqyechim_array)
print("--------------------------------------------")
print("taqribiy yechimlar")
print(find_U_x_l_taqribiy_yechim_array)
print("---------------------------------------------")
print("xatolik")
for i in range(0, N + 1):
    print("[", find_aniqyechim_array[i], " - ", find_U_x_l_taqribiy_yechim_array[i], "] = ", find_xatoliklar_array[i])
plt.plot(find_X_l_array,find_aniqyechim_array, 'o-',  label = "Aniq yechim", color = "green")
plt.plot(find_X_l_array,find_U_x_l_taqribiy_yechim_array,'o--', label = "Taqribiy yechim", color = "red")
plt.title("Aniq va Taqribiy yechimlar grafigi")
plt.xlabel(f"N = {N},  A = {A}", loc = "right", color = "red", size = 15)
plt.ylabel("natijalar")
plt.legend(loc = "best", fontsize = 12)
plt.grid(color = "blue", linestyle = ":", linewidth = 0.5)
plt.show()