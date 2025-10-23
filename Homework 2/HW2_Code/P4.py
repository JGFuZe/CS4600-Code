# CS4600 - Homework 2 - Problem 4
# Jonah Gallagher

import math
import numpy as np
import matplotlib.pyplot as plt

# Defines the cost function f(x, y) = sin(x) + cos(y) + ((x - y)^2)/4
def inefficiencyCost(x, y):
    return ((math.sin(x)) + (math.cos(y)) + (((x - y)**2) / 4))

# inefficiencyCost End

# Gradient components computed by hand
def inefficiencyGrad(x, y):
    # Gradient x-component = cos(x) + 0.5*(x - y)
    gradX = ((math.cos(x)) + ((x - y) / 2))

    # Gradient y-component = -sin(y) - 0.5*(x - y)
    gradY = (((-1) * math.sin(y)) - ((x - y) / 2))

    # Return the gradient vector as a tuple
    return (gradX, gradY)

# Compute the Euclidean norm of the gradient vector
def gradientNorm(gradVector):

    # Use sqrt(gx^2 + gy^2) to measure gradient magnitude
    return (math.sqrt(((gradVector[0]**2) + (gradVector[1]**2))))
# gradientNorm End

# Perform the steepest descent iteration with a fixed step size
def steepestDescent(startPoint, stepSize, tolerance, maxIterations):

    # Split out the starting coordinates
    xCurrent, yCurrent = startPoint

    # Track every point for plotting later
    pathPoints = [(xCurrent, yCurrent)]

    # Grab the gradient at the starting location
    gradVector = inefficiencyGrad(xCurrent, yCurrent)

    # Measure how large the gradient is
    gradNorm = gradientNorm(gradVector)

    # Start the iteration counter at zero
    iterations = 0

    # Keep iterating while we have a big gradient and room for more iterations
    while (((gradNorm > tolerance)) and (iterations < maxIterations)):
        # Pull apart the gradient components
        gradX, gradY = gradVector

        # Newton style update using the fixed step size alpha
        xNext = (xCurrent - (stepSize * gradX))

        # Do the same for the y value
        yNext = (yCurrent - (stepSize * gradY))

        # Move the current point to the newly computed point
        xCurrent = xNext
        yCurrent = yNext

        # Count that we just finished an iteration
        iterations += 1

        # Re-evaluate the gradient at the updated point
        gradVector = inefficiencyGrad(xCurrent, yCurrent)

        # Update the gradient size for the stopping rule
        gradNorm = gradientNorm(gradVector)

        # Store this point so our plot shows the full path
        pathPoints.append((xCurrent, yCurrent))

    # Get the final function value once we stop iterating
    finalValue = inefficiencyCost(xCurrent, yCurrent)

    # Return Array of values [iterations, xFinal, yFinal, finalFunctionValue, finalGradientNorm, pathPoints]
    return ([iterations, xCurrent, yCurrent, finalValue, gradNorm, pathPoints])

# steepestDescent End

# Plot the contour of the function along with the optimization path
def plotDescentPath(pathPoints, bounds):
    # Unpack the lower and upper values for both axes
    lower, upper = bounds

    # Grab just the x coordinates from the stored path
    pathX = [point[0] for point in pathPoints]

    # Grab just the y coordinates from the stored path
    pathY = [point[1] for point in pathPoints]

    # Compute the extreme values reached by the path so we can expand the view
    xMin = min(pathX)
    xMax = max(pathX)
    yMin = min(pathY)
    yMax = max(pathY)

    # Add a little padding so the path is not plotted right on the border
    padding = 0.5
    xLower = min(lower, (xMin - padding))
    xUpper = max(upper, (xMax + padding))
    yLower = min(lower, (yMin - padding))
    yUpper = max(upper, (yMax + padding))

    # Build the x values across the adjusted plotting range
    gridX = np.linspace(xLower, xUpper, 200)

    # Build the y values across the adjusted plotting range
    gridY = np.linspace(yLower, yUpper, 200)

    # Create the mesh grid needed for contour plots
    meshX, meshY = np.meshgrid(gridX, gridY)

    # Evaluate f(x, y) over the whole mesh grid
    meshZ = ((np.sin(meshX)) + (np.cos(meshY)) + (((meshX - meshY)**2) / 4))

    # Start a new figure window
    plt.figure(figsize=(8, 6))

    # Pick how many contour levels we want
    contourLevels = 25

    # Draw the contour plot using the viridis color map
    plt.contour(meshX, meshY, meshZ, contourLevels, cmap="viridis")

    # Plot the path with red circles connecting each step
    plt.plot(pathX, pathY, marker="o", color="red", linewidth=1.5, markersize=4, label="Path")

    # Mark the final point with a black X for clarity
    plt.scatter(pathX[-1], pathY[-1], color="black", marker="x", s=50, label="Final Point")

    # Lock the axis limits so the enlarged window stays in view
    plt.xlim(xLower, xUpper)
    plt.ylim(yLower, yUpper)

    # Label the plot elements for the assignment write up
    plt.title("Steepest Descent Path for f(x, y)")
    plt.xlabel("x")
    plt.ylabel("y")

    # Show a legend so the markers are easy to see
    plt.legend()

    # Add a light grid so the path is easy to see
    plt.grid(True, linestyle="--", alpha=0.3)

    # Show the plot
    plt.show()

    # Close the figure so we do not stack plots between runs
    plt.close()

# plotDescentPath End

# Main driver to execute the required parts of Problem 4
def main():
    
    # Given step size in assignment
    alpha = 0.2

    # given Starting point in assignment
    startPoint = (1.0, 2.0)

    # Part (c): Full steepest descent run with tolerance and max iterations
    tolerance = (10**(-4))
    # Cap the run so we do not loop forever
    maxIterations = 10000

    # Call steepest descent and capture all returned values
    descentResult = steepestDescent(startPoint, alpha, tolerance, maxIterations)

    # Grab each value out of the returned array
    descentIterations = descentResult[0]
    descentX = descentResult[1]
    descentY = descentResult[2]
    descentValue = descentResult[3]
    descentGrad = descentResult[4]
    descentPath = descentResult[5]

    # Print the summary information for the full run
    print("Steepest Descent Summary:")
    print(f"Iterations: {descentIterations}")
    print(f"Final Point: ({descentX:.6f}, {descentY:.6f})")
    print(f"Final Function Value: {descentValue:.6f}")
    print(f"Final Gradient Norm: {descentGrad:.6f}\n")

    # Generate the contour plot with the optimization path overlay
    plotBounds = (-3, 3)
  
    # Build the plot and show it
    plotDescentPath(descentPath, plotBounds)

# Main End

if (__name__ == "__main__"):
    main()
