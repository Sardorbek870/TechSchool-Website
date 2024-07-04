
import numpy as np
pi = np.pi
e = np.e
# pi = 3.14
# e = 2.71

N = int(input("N = "))
M = int(input("M = "))
A = float(input("A = "))

l_values = np.arange(N+1)
x_l_array = 1 * np.cos(pi * l_values / N)

print("x(l) array ")
print(x_l_array)
print("-" * 30)

k_values = np.arange(M+1)
y_k_array = 1 * np.cos(pi * k_values / M)

print("y(k) array")
print(y_k_array)
print("-" * 30)


T_n_x_l_matrix = np.zeros((N+1, N+1))
for n in range(N+1):
    for l in range(N+1):
        T_n_x_l_matrix[n,l] = np.cos((n * pi * l) / N)

print("T_n(x_l) matritsaning qiymatlari")
print(T_n_x_l_matrix)
print("-" * 30)

T_n_y_k_matrix = np.zeros((M+1, M+1))
for n in range(N+1):
    for k in range(M+1):
        T_n_y_k_matrix[n,k] = np.cos((n * pi * k) / M)

print("T_n(y_k) matritsaning qiymatlari")
print(T_n_y_k_matrix)
print("-" * 30)

aniq_yechim_U_x_l = np.zeros((N+1, M+1))
for i in range(N+1):
    for j in range(M+1):
        aniq_yechim_U_x_l[i,j] = (1 - pow(x_l_array[i], 2)) * (1 - pow(y_k_array[j], 2)) * pow(e, (A * x_l_array[i] * y_k_array[j]))
aniq_yechim_U_x_l = aniq_yechim_U_x_l.flatten()
print("U(x,y) aniq yechim")
print(aniq_yechim_U_x_l)
print("-" * 30)


def find_f_x_y_matrix(x,y):
    return (2 * np.exp(A * x * y) * (2 - pow(x,2) - pow(y,2)) *
            (1 + 2 * A * x * y) - (1 - pow(x,2)) * (1 - pow(y,2)) *
            pow(A,2) * np.exp(A * x * y) * (pow(x,2) + pow(y,2)))


f_x_y_k_matrix = np.zeros((N+1, M+1))
for i in range(N+1):
    for j in range(M+1):
        f_x_y_k_matrix[i,j] = find_f_x_y_matrix(x_l_array[i], y_k_array[j])

print("f(x(l), y(k)) matrix")
print(f_x_y_k_matrix)
print("-" * 30)


def find_g_i_j(i,j):
    C_i = 1
    if i == 0 or i == N:
        C_i = 2

    C_j = 1
    if j == 0 or j == M:
        C_j = 2
    result = 4 / (M * N * C_i * C_j)
    sum = 0

    for l in range(N+1):
        for k in range(M+1):
            C_l = 1
            if(l == 0 or l == N):
                C_l = 2
            C_k = 1
            if k == 0 or k == M:
                C_k = 2

            sum += 1 / (C_l * C_k) * f_x_y_k_matrix[l, k] * T_n_x_l_matrix[i,l] * T_n_y_k_matrix[j,k]
    result *= sum
    return result


g_i_j_matrix = np.zeros((N+1, M+1))
for i in range(N+1):
    for j in range(M+1):
        g_i_j_matrix[i,j] = find_g_i_j(i,j)

print("g_i_j_matrix")
print(g_i_j_matrix)
print("-" * 30)

# 4 a)
A_matrix = np.zeros(((M+1)*(N+1), (M+1)*(N+1)))
b_array = np.zeros((N+1)*(M+1))


index_of_n = 0

for i in range(2, N+1):
    for j in range(2, M+1):
        A_temporary_matrix = np.zeros((N + 1, M + 1))
        A_temporary_matrix[i, j - 2] += 4 * i * (i**2 - 1) * (j + 1)
        A_temporary_matrix[i, j] += 4 * i * (i**2 - 1) * (-2 * j)

        if j != M and j != (M - 1):
            A_temporary_matrix[i, j + 2] += 4 * i * (pow(i,2) - 1) * (j - 1)

        A_temporary_matrix[i - 2, j] += 4 * j * (pow(j,2) - 1) * (i+1)
        A_temporary_matrix[i,j] += 4 * j * (pow(j,2) - 1) * (-2 * i)

        if i != N and i != N - 1:
            A_temporary_matrix[i + 2, j] += (4 * j * (pow(j,2) - 1) * (i - 1))

        A_matrix[index_of_n, :] = A_temporary_matrix.flatten()


# b_array[index_of_n] = -1 * ((i + 1) * ((j + 1) * g_i_j_matrix[i - 2, j - 2] - 2 * j * g_i_j_matrix[i - 2, j] +
        #                                        (j - 1) * g_i_j_matrix[i - 2, j + 2]) -
        #                                        2 * i * ((j + 1) * g_i_j_matrix[i, j - 2] - 2 * j * g_i_j_matrix[i,j]  +
        #                                                 (j - 1) * g_i_j_matrix[i, j + 2]) +
        #                                                  (i - 1) * ((j + 1) * g_i_j_matrix[i+2, j-2] -
        #                                                 2 * j * g_i_j_matrix[i+2, j] +
        #                                                 (j - 1) * g_i_j_matrix[i+2, j+2]))
        if i+2 >= N and j+2 >=M:
            b_array[index_of_n] = -1 * ((i+1) * ((j + 1) * g_i_j_matrix[i-2, j-2] - 2 * j * g_i_j_matrix[i-2, j])-
                                    2 * i * ( (j + 1) * g_i_j_matrix[i, j-2] - 2*j*g_i_j_matrix[i,j]))
        elif i+2 >= N:
            b_array[index_of_n] = -1 * ((i+1) * ((j+1) * g_i_j_matrix[i-2, j - 2] - 2 * j * g_i_j_matrix[i-2, j] + (j-1) * g_i_j_matrix[i-2, j+2]) -
                  2 * i * ((j + 1) * g_i_j_matrix[i,j-2] - 2 * j * g_i_j_matrix[i,j] + (j - 1) * g_i_j_matrix[i, j+2]))
        elif j+2 >= M:
            b_array[index_of_n] = -1 * ((i+1) * ((j + 1) * g_i_j_matrix[i-2, j-2] - 2 * j * g_i_j_matrix[i-2, j])-
                  2 * i * ((j+1) * g_i_j_matrix[i,j-2] - 2 * j * g_i_j_matrix[i,j] )+
                  (i-1) * ((j+1) * g_i_j_matrix[i+2, j-2] - 2 * j * g_i_j_matrix[i+2,j]))
        else:
            b_array[index_of_n] = -1 * ((i+1) * ((j+1) * g_i_j_matrix[i-2, j-2] - 2 * j * g_i_j_matrix[i-2, j] + (j-1) * g_i_j_matrix[i-2, j-2]) -
                  2 * i * ((j+1) * g_i_j_matrix[i, j - 2] - 2 * j * g_i_j_matrix[i,j] + (j - 1) * g_i_j_matrix[i,j+2]) +
                  (i - 1) * ((j+1) * g_i_j_matrix[i+2, j-2] - 2 * j * g_i_j_matrix[i+2, j] + (j-1) * g_i_j_matrix[i+2, j+2]))



        index_of_n += 1


#4 b)
# index_of_n = (N - 1) * (M - 1)
for j in range(2, M+1):
    A_temporary_matrix = np.zeros((N + 1, M + 1))
    for i in range(1, N+1):
        if i % 2 == 1:
            A_temporary_matrix[i,j] = 1

    A_matrix[index_of_n, :] = A_temporary_matrix.flatten()
    index_of_n += 1


#4 c)

for j in range(2, M+1):
    A_temporary_matrix = np.zeros((N + 1, M + 1))
    for i in range(0, N+1):
        if i % 2 == 0:
            if i == 0 and j == 0:
                A_temporary_matrix[i, j] = 1 / 4
            elif i == 0:
                A_temporary_matrix[i, j] = 1 / 2
            else:
                A_temporary_matrix[i, j] = 1

    A_matrix[index_of_n, :] = A_temporary_matrix.flatten()
    index_of_n += 1

#4 d)



for i in range(0, N+1):
    A_temporary_matrix = np.zeros((N + 1, M + 1))
    for j in range(1, M+1):
        if j % 2 == 1:
            A_temporary_matrix[i,j] = 1

    A_matrix[index_of_n, :] = A_temporary_matrix.flatten()
    index_of_n += 1

#4 e)

for i in range(0, N+1):
    A_temporary_matrix = np.zeros((N + 1, M + 1))
    for j in range(0, M+1):
        if i == 0 and j == 0:
            A_temporary_matrix[i,j] = 1 / 4
        elif j == 0:
            A_temporary_matrix[i, j] = 1 / 2
        elif j % 2 == 0:
            A_temporary_matrix[i, j] = 1
    A_matrix[index_of_n, :] = A_temporary_matrix.flatten()
    index_of_n += 1

print("A_matrix")
print(A_matrix.shape)
print(A_matrix)
print("-------------------------")
filename = "file.txt"
np.savetxt(filename, A_matrix)

print(f"b array in A * x = b  size of b = {b_array.shape}")
print(b_array)
print("-------------------------")

gauss_result = np.linalg.solve(A_matrix, b_array)
print("gauss result")
print(gauss_result)
print(gauss_result.size)
print("-" * 30)


taqribiy_yechim_U_x_y = np.zeros(((N+1) * (M+1)))
index = 0


for l in range(len(l_values)):
    for k in range(len(k_values)):
        result = 0
        index_to_gauss = 0
        for i in range(N+1):
            for j in range(M+1):
                multipliedValue = 1
                if i == 0 and j == 0:
                    multipliedValue = 1 / 4
                elif i == 0:
                    multipliedValue = 1 / 2
                elif j == 0:
                    multipliedValue = 1 / 2

                #result += gauss_result[i * N + j] * T_n_x_l_matrix[i,l] * T_n_y_k_matrix[j,k]
                result += multipliedValue * gauss_result[index_to_gauss] * T_n_x_l_matrix[i, l] * T_n_y_k_matrix[j, k]
                index_to_gauss += 1

        taqribiy_yechim_U_x_y[index] = result
        index = index + 1
print("taqribiy yechim U(x,y)")
print(taqribiy_yechim_U_x_y)
print("---------------------")

hatoliklar_array = np.zeros(len(aniq_yechim_U_x_l))
for i in range(len(hatoliklar_array)):
    hatoliklar_array[i] = np.fabs(aniq_yechim_U_x_l[i] - taqribiy_yechim_U_x_y[i])

print("hatolik")
for i in range(len(aniq_yechim_U_x_l)):
    print(f"[{aniq_yechim_U_x_l[i]}  -  {taqribiy_yechim_U_x_y[i]}  =  {hatoliklar_array[i]}")