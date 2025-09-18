import numpy as np
import matplotlib.pyplot as plt

# Define the function
def f(x):
    return (x**3) - 6*(x**2) + 11*x - 6.1

# Generate x values
x = np.linspace(-10, 10, 100)

# Compute y values
y = f(x)

# Plot the function
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of f(x) = x^3 - 6x^2 + 11x - 6.1')
plt.grid(True)
plt.show()