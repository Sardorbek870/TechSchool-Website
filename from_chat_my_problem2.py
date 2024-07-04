import numpy as np
import matplotlib.pyplot as plt

# Define parameters
N = int(input("N = "))
eps = 0.01  # Given in the problem

# Chebyshev points
l_values = np.arange(N + 1)
y_l_array = -np.cos(np.pi * l_values / N)

# Function to evaluate the right-hand side of the differential equation
def find_f(y_l):
    return (1 / 8) * (y_l + 1)

# Create the matrix of Chebyshev polynomials evaluated at the points
find_T_n_matrix = np.zeros((N + 1, N + 1))
find_T_n_matrix[0, :] = 1
find_T_n_matrix[1, :] = y_l_array
for n in range(2, N + 1):
    find_T_n_matrix[n, :] = 2 * y_l_array * find_T_n_matrix[n - 1, :] - find_T_n_matrix[n - 2, :]

# Right-hand side evaluated at the Chebyshev points
find_f_massiv = find_f(y_l_array)

# Construct the b vector
def find_b_n(n):
    c_n = 2 if (n == 0 or n == N) else 1
    sum = 0
    for i in range(N + 1):
        c_l = 2 if (i == 0 or i == N) else 1
        sum += (1 / c_l) * find_f_massiv[i] * find_T_n_matrix[n, i]
    return 2 / (N * c_n) * sum

free_b_array = np.array([find_b_n(i) for i in range(N + 1)])

# Construct the A matrix
a_matrix = np.zeros((N + 1, N + 1))
for i in range(N + 1):
    if i == 0:
        c_i = 2
    else:
        c_i = 1

    for p in range(i + 2, N + 1):
        if (p - i) % 2 == 0:
            a_matrix[i, p] += (1 / c_i) * p * (p**2 - i**2)
    a_matrix[i, i] -= 1 / (2 * c_i)

# Add boundary conditions to the matrix
for i in range(N + 1):
    a_matrix[N - 1, i] = 1  # U(-1) = 0
    a_matrix[N, i] = 1 if i % 2 == 0 else -1  # U(1) = 0

# Solve the system
a_gauss_result = np.linalg.solve(a_matrix, free_b_array)

# Define the approximate solution
def find_U_taqribiy(y_l):
    sum = 0
    for i in range(N + 1):
        sum += a_gauss_result[i] * find_T_n_matrix[i, y_l]
    return sum

find_U_taqribiy_array = np.array([find_U_taqribiy(i) for i in range(N + 1)])

# Define the exact solution
def find_U_aniq_yechim(y):
    term1 = (eps - 0.5) / (1 - np.exp(-1 / eps)) * (1 - np.exp(-(y + 1) / (2 * eps)))
    term2 = -eps * (y + 1) / 2
    term3 = (y + 1)**2 / 8
    return term1 + term2 + term3

find_U_aniq_yechim_massiv = np.array([find_U_aniq_yechim(y) for y in y_l_array])

# Calculate the error
find_xatolik_massiv = np.abs(find_U_taqribiy_array - find_U_aniq_yechim_massiv)

# Print results
print("y(l) larning qiymati")
for i in range(N + 1):
    print(f"y({i}) = {y_l_array[i]}")
print("-------------------------------------------------")

print("T larning qiymatlari")
for i in range(N + 1):
    for j in range(N + 1):
        print(f"T{i}(y({j})) = {find_T_n_matrix[i, j]}")
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
for i in range(N + 1):
    print(f"U(y({i})) = {find_U_taqribiy_array[i]}")
print("-" * 20)

print("Aniq yechimlar")
for i in range(N + 1):
    print(f"U(y({i})) = {find_U_aniq_yechim_massiv[i]}")
print("-" * 20)

print("xatoliklar")
for i in range(N + 1):
    print(f"{find_U_taqribiy_array[i]} - {find_U_aniq_yechim_massiv[i]} = {find_xatolik_massiv[i]}")
print("-" * 20)

print(f"xatoliklarning o'rta arifmetigi = {np.mean(find_xatolik_massiv)}")

# Plot results
plt.plot(y_l_array, find_U_aniq_yechim_massiv, 'o--', label="Aniq yechim", color="green")
plt.plot(y_l_array, find_U_taqribiy_array, 'o--', label="Taqribiy yechim", color="red")
plt.title("Aniq va Taqribiy yechimlar grafigi")
plt.xlabel(f"N = {N}, eps = {eps}", loc="right", color="red", size=15)
plt.ylabel("natijalar")
plt.legend(loc="best", fontsize=12)
plt.grid(color="blue", linestyle=":", linewidth=0.5)
plt.show()
