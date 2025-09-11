import numpy as np
import matplotlib.pyplot as plt

# Define the function
def f(x):
    return ((-12) - (21*x) + (18)*(x**2) - (2.75)*(x**3))

# Generate x values
x = np.linspace(-10, 10, 100)

# Compute y values
y = f(x)

# Plot the function
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of f(x) = -12 - 21x + 18x^2 - 2.75x^3')
plt.grid(True)
plt.show()

