import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
e = np.e
N = int(input("N = "))
eps = float(input("eps = "))

l_values = np.arange(N+1)
y_l_array =-1 * np.cos(pi * l_values / N)

find_T_n_matrix = np.zeros((N + 1, N + 1))
find_T_n_matrix[0, :] = 1
find_T_n_matrix[1, :] = y_l_array
for n in range(2, N + 1):
    find_T_n_matrix[n, :] = 2 * y_l_array * find_T_n_matrix[n - 1, :] - find_T_n_matrix[n - 2, :]


def find_f(y_l):
    return (1 / 8) * (y_l + 1)

find_f_massiv = np.zeros(N+1)

for i in range(0,N+1):
    find_f_massiv[i] = find_f(y_l_array[i])

def find_b_n(n):
    if (n == 0 or n == N):
        c_n = 2
    else:
        c_n = 1

    sum = 0
    for i in range(0, N + 1):
        if (i == 0 or i == N) :
            c_l = 2
        else:
            c_l = 1
        sum += (1 / c_l) * find_f_massiv[i] * find_T_n_matrix[n, i]

    return 2 / (N * c_n) * sum

a_matrix = np.zeros((N + 1, N + 1))
free_b_array = np.zeros(N + 1)

for i in range(0, N - 1):
    if i == 0 or i == N:
        c_i = 2
    else:
        c_i = 1

    for p in range(i+2, N+1):
        if (p - i) % 2 == 0:
            a_matrix[i,p] += eps / c_i * p * (pow(p,2) - pow(i, 2))
    a_matrix[i,i] -= 1/(2)
    free_b_array[i] = find_b_n(i)


for i in range(0, N+1):
    a_matrix[N,i] = 1

    if i % 2 == 0:
        a_matrix[N-1, i] = 1
    else:
        a_matrix[N-1, i] = -1


a_gauss_result = np.linalg.solve(a_matrix, free_b_array)

def find_U_taqribiy(y_l):
    sum = 0
    for i in range(0, N+1):
        sum += a_gauss_result[i] * find_T_n_matrix[i, y_l]
    return sum

find_U_taqribiy_array = np.zeros(N+1)

for i in range(N+1):
    find_U_taqribiy_array[i] = find_U_taqribiy(i)



def find_U_aniq_yechim(y_l):

    return ((eps - 0.5) / (1 - pow(e, ((-1 / eps)))) *
            (1 - pow(e, (-(y_l_array[y_l] + 1) / (2 * eps)))) - eps * ((y_l_array[y_l] + 1) / 2) +
            (pow((y_l_array[y_l] + 1), 2) / 8))



find_U_aniq_yechim_massiv = np.zeros(N + 1)
for i in range(N + 1):
    find_U_aniq_yechim_massiv[i] = find_U_aniq_yechim(i)

find_xatolik_massiv = np.zeros(N + 1)

for i in range(N + 1):
    find_xatolik_massiv[i] = np.fabs(find_U_aniq_yechim_massiv[i] - find_U_taqribiy_array[i])

print("y(l) larning qiymati")
for i in range(0, N + 1):
    print("y(", i, ")=", y_l_array[i])
print("-------------------------------------------------")

print("T larning qiymatlari")
for i in range(0, N + 1):
    for j in range(0, N + 1):
        print("T", i, "(y(", j, ")) = ", find_T_n_matrix[i, j])
    print("")
print("---------------------------------------------------")

print("Ax = b dagi A matritsa")
print(a_matrix)
print("-" * 20)

print("Ax = b dagi b vektor ning qiymatlari")
print(free_b_array)
print("-" * 20)

print("Ax = b tenglamani gauss usuli orqali hisoblangan yechimlari")
print(a_gauss_result)
print("-" * 20)

print("Taqribiy yechimlar")
for i in range(N+1):
    print(f"U(y({i})) = {find_U_taqribiy_array[i]}")
print("-" * 20)

print("Aniq yechimlar")
for i in range(N+1):
    print(f"U(y({i})) = {find_U_aniq_yechim_massiv[i]}")
print("-" * 20)

print("xatoliklar")
for i in range(N+1):
    print(f"{find_U_taqribiy_array[i]} - {find_U_aniq_yechim_massiv[i]} = {find_xatolik_massiv[i]}")
print("-" * 20)

print(f"xatoliklarning o'rta arifmetigi = {sum(find_xatolik_massiv) / N}")

plt.plot(y_l_array, find_U_aniq_yechim_massiv, 'o--', label="Aniq yechim", color="green")
plt.plot(y_l_array, find_U_taqribiy_array, 'o--', label="Taqribiy yechim", color="red")
plt.title("Aniq va Taqribiy yechimlar grafigi")
plt.xlabel(f"N = {N},  eps = {eps}", loc="right", color="red", size=15)
plt.ylabel("natijalar")
plt.legend(loc="best", fontsize=12)
plt.grid(color="blue", linestyle=":", linewidth=0.5)

plt.show()

c=list([1 for i in range(N+1)])
c[0]=c[N]=2
# Birinchi tartibli Hosilani hisoblash
a_1_array = np.zeros(N+1)
for n in range(N+1):
    for p in range(n+1, N+1):
        if (p+n - 1) % 2 == 0:
            a_1_array[n] += 2 * p * a_gauss_result[p] / c[n]

U_y_hosila_1 = np.zeros(N+1)
for l in range(N+1):
    for i in range(N+1):
        U_y_hosila_1[l] += a_1_array[i] * find_T_n_matrix[i,l]

U_y_hosila_1_aniq_yechim = np.zeros(N+1)
for i in range(N+1):
    U_y_hosila_1_aniq_yechim[i] = (eps - 1/2) * np.exp(-1 * (y_l_array[i] + 1) / (2 * eps)) / (2 * eps * (1 - np.exp(-1 / eps))) + (y_l_array[i] + 1) / 4 - eps / 2

# Ikkinchi tartibli hosila
a_2_array = np.zeros(N+1)
for n in range(N+1):
    for p in range(n+2, N+1):
        if (p - n) % 2 == 0:
            a_2_array[n] += p * (pow(p,2) - pow(n,2)) * a_gauss_result[p] / c[n]

U_y_hosila_2 = np.zeros(N+1)
for l in range(N+1):
    for i in range(N+1):
        U_y_hosila_2[l] += a_2_array[i] * find_T_n_matrix[i,l]

U_y_hosila_2_aniq_yechim = np.zeros(N+1)
for i in range(N+1):
    U_y_hosila_2_aniq_yechim[i] = 1/4 - (eps - 1/2) * np.exp((y_l_array[i] + 1) / (-2 * eps)) / (4 * pow(eps,2) * (1 - np.exp(-1 / eps)))

print("1-tartibli hosila taqribiy - aniq = xatolik")

for i in range(N+1):
    print(f"{U_y_hosila_1[i]} - {U_y_hosila_1_aniq_yechim[i]} = {np.fabs(U_y_hosila_1[i] - U_y_hosila_1_aniq_yechim[i])}")
print("-" * 20)


plt.plot(U_y_hosila_1,"o-", label = "Taqribiy yechim", color = "red")
plt.plot(U_y_hosila_1_aniq_yechim, "o-", label = "Aniq yechim", color = "blue")
plt.title("Birinchi tartibli hosila")
plt.ylabel('Hosila')
plt.xlabel('l ning qiymatlari')
plt.legend(loc='best', fontsize=10)
plt.show()


print("2-tartibli hosila taqribiy - aniq = xatolik")

for i in range(N+1):
    print(f"{U_y_hosila_2[i]} - {U_y_hosila_2_aniq_yechim[i]} = {np.fabs(U_y_hosila_2[i] - U_y_hosila_2_aniq_yechim[i])}")
print("-" * 20)

plt.plot(U_y_hosila_2,"-o", label = "Taqribiy yechim", color = "red")
plt.plot(U_y_hosila_2_aniq_yechim, "-o", label = "Aniq yechim", color = "blue")
plt.title("Ikkinchi tartibli hosila")
plt.ylabel("Holisa")
plt.xlabel('l ning qiymatlari')
plt.legend(loc='best', fontsize=10)
plt.show()

