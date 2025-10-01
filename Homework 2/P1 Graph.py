# CS4600 - Homework 2 - Problem 1 Graph
# Jonah Gallagher

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define the function
def f(x, y):
    # Return f(x, y) = x^3 + y^3 - 3xy
    return (x**3 + y**3 - (3 * x * y))


# 
def main():
    # Init Variables
    X, Y, Z = None, None, None
    fig = None
    axis = None
    surf = None

    # Set Grid points
    gridPoints = 200

    # Set boundries for x and y
    x = np.linspace(-2, 2, gridPoints)
    y = np.linspace(-2, 2, gridPoints)

    # Create meshgrid
    X, Y = np.meshgrid(x, y)

    # Compute Z values
    Z = f(X, Y)

    # Plot the surface
    fig = plt.figure(figsize=(8, 6))

    # 3D axis
    axis = fig.add_subplot(111, projection="3d")
 
    # Create the surface plot 
    surf = axis.plot_surface(X, Y, Z, cmap=cm.viridis, edgecolor="none")

    # Labels and title
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_zlabel("f(x, y)")
    axis.set_title("Surface plot of f(x, y) = x^3 + y^3 - 3xy")

    # Color bar
    fig.colorbar(surf, shrink=0.6, aspect=12, label="f(x, y)")

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
