import numpy as np
import matplotlib.pyplot as plt

# Parameters
eps = 0.01  # Given in the problem
N = 50  # Number of grid points

# Chebyshev nodes and differentiation matrix
l_values = np.arange(N + 1)
y_l_array = -np.cos(np.pi * l_values / N)

# Chebyshev differentiation matrix
D = np.zeros((N + 1, N + 1))
for i in range(N + 1):
    for j in range(N + 1):
        if i != j:
            D[i, j] = (-1)**(i + j) / (y_l_array[i] - y_l_array[j])
        elif i != 0:
            D[i, i] = -y_l_array[i] / (2 * (1 - y_l_array[i]**2 + 1e-12))
        else:
            D[i, i] = (2 * (i**2) + 1) / 6

# Differential operators in the spectral method
A = eps * D.dot(D) - 0.5 * D

# Right-hand side vector
b = 1 / 8 * (y_l_array + 1)

# Boundary conditions
A[0, :] = 0
A[0, 0] = 1
A[N, :] = 0
A[N, N] = 1

# Solve for coefficients
U_coeffs = np.linalg.solve(A, b)

# Function to evaluate approximate solution at y
def U_approx(y):
    return sum(U_coeffs[n] * np.cos(n * np.arccos(y)) for n in range(N + 1))

# Exact solution function
def U_exact(y):
    return ((eps - 0.5) / (1 - np.exp(-1 / eps))) * (1 - np.exp(-(y + 1) / (2 * eps))) - eps * ((y + 1) / 2) + ((y + 1)**2 / 8)

# Evaluate solutions over a range of y values
y_values = np.linspace(-1, 1, 200)
U_approx_values = np.array([U_approx(y) for y in y_values])
U_exact_values = np.array([U_exact(y) for y in y_values])

# Plotting
plt.plot(y_values, U_approx_values, label='Approximate Solution (Spectral Method)', color='blue')
plt.plot(y_values, U_exact_values, label='Exact Solution', linestyle='--', color='red')
plt.title('Comparison of Spectral Method Approximation with Exact Solution')
plt.xlabel('y')
plt.ylabel('U(y)')
plt.legend()
plt.grid(True)
plt.show()

# Error calculation
error = np.linalg.norm(U_approx_values - U_exact_values, np.inf)
print(f"Maximum error: {error}")

