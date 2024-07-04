import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 10  # Number of grid points
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

# Construct the A matrix and b vector
A = np.zeros((N + 1, N + 1))
b = np.zeros(N + 1)

# Filling the A matrix and b vector according to the differential equation
for i in range(1, N):
    for j in range(N + 1):
        if j == 0 or j == N:
            c_j = 2
        else:
            c_j = 1
        A[i, j] = (eps * (j**2) * (1 if i == j else 2 * (-1)**(i-j))) / c_j - 0.5 * find_T_n_matrix[j, i]
    b[i] = find_f_massiv[i]

# Boundary conditions
A[0, 0] = A[0, N] = 1
A[N, 0] = A[N, N] = (-1)**(N)

# Solve the system with regularization
reg_param = 1e-6 * np.eye(N + 1)  # Small regularization parameter
a = np.linalg.solve(A + reg_param, b)

# Define the approximate solution
def find_U_taqribiy(y):
    return sum(a[j] * np.cos(j * np.arccos(y)) for j in range(N + 1))

find_U_taqribiy_array = np.array([find_U_taqribiy(y) for y in y_l_array])

# Define the exact solution
def find_U_aniq_yechim(y):
    return ((eps - 0.5) / (1 - np.exp(-1 / eps))) * (1 - np.exp(-(y + 1) / (2 * eps))) - eps * (y + 1) / 2 + (y + 1)**2 / 8

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
print(A)
print("-" * 20)

print("Ax = b dagi b vektor ning qiymatlari")
print(b)
print("-" * 20)

print("Ax = b tenglamani gauss usuli orqali hisoblangan yechimlari")
print(a)
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
