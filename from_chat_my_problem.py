import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve


def chebyshev_points(N):
    return np.cos(np.pi * np.arange(N + 1) / N)


def chebyshev_differentiation_matrix(N):
    x = chebyshev_points(N)
    c = np.ones((N + 1,))
    c[0] = 2
    c[-1] = 2
    c = c * (-1) ** np.arange(N + 1)

    X = np.tile(x, (N + 1, 1))
    dX = X - X.T

    D = (c[:, None] / c[None, :]) / (dX + np.eye(N + 1))
    D = D - np.diag(np.sum(D, axis=1))

    return D


def solve_diff_eq(N):
    x = chebyshev_points(N)
    D = chebyshev_differentiation_matrix(N)
    D2 = D @ D

    # Right-hand side of the equation
    f = (1 / 8) * (x + 1)

    # Matrix A representing the equation
    A = 0.01 * D2 - 0.5 * np.eye(N + 1)

    # Apply boundary conditions: U(-1) = 0, U(1) = 0
    A[0, :] = 0
    A[-1, :] = 0
    A[0, 0] = 1
    A[-1, -1] = 1
    f[0] = 0
    f[-1] = 0

    # Solve the linear system Au = f
    U = solve(A, f)

    return x, U


def exact_solution(y):
    C2 = 1 / (2 * (np.exp(-np.sqrt(50)) - np.exp(3 * np.sqrt(50))))
    C1 = -C2 * np.exp(2 * np.sqrt(50))
    return C1 * np.exp(np.sqrt(50) * y) + C2 * np.exp(-np.sqrt(50) * y) - (1 / 4) * y - (1 / 4)


# Parameters
N = 32  # Number of Chebyshev points

# Solve the differential equation
x, U_numerical = solve_diff_eq(N)

# Compute the exact solution at the Chebyshev points
U_exact = exact_solution(x)

# Plot the solutions
plt.plot(x, U_numerical, 'o-', label='Numerical solution')
plt.plot(x, U_exact, 'x-', label='Exact solution')
plt.xlabel('y')
plt.ylabel('U(y)')
plt.title('Comparison of Numerical and Exact Solutions')
plt.legend()
plt.grid(True)
plt.show()

# Compute the error
error = np.abs(U_numerical - U_exact)
print("Error at Chebyshev points:", error)
