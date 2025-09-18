# Defines the function for root finding
# -X^2 - 2
def rootFunc(x):
    return ((1)*(x**2) - 2)


def checkRootExists(lBound, rBound):
    # value to return
    returnVal = 0

    #
    rootLower = (rootFunc(lBound))

    # if f(x) Left Bound Value * Right Bound Value is < 0 Then There is a root
    if ((rootLower * rootFunc(rBound)) < 0):
        returnVal = -1

    elif ((rootLower) == 0): # Root is the lower bound
        returnVal = 0

    else: # If not then there is no garenteed root
        returnVal = 1

    return returnVal

def bisectionMethod(lBound, rBound, tolerance, maxIterations, rootFunc):
    # Initialize Variables
    midpoint = 0
    previousMid = 0
    value = 0

    # Set initial bounds
    xLeft = lBound
    xRight = rBound

    # iterations counter
    iterations = 0
    
    # Check Error
    relError = abs(xRight - xLeft)

    # check tolerances
    while ((relError > tolerance) and (iterations < maxIterations)):
        # Increment Iterations
        iterations += 1

        # Get Midpoint
        midpoint = ((xLeft + xRight) / 2) 

        # Calculate value at midpoint
        value = rootFunc(midpoint)

        if ((rootFunc(xLeft) * value) < 0):
            xRight = midpoint
        elif ((rootFunc(xLeft) * value) > 0):
            xLeft = midpoint
        else:
            # Root is found
            break

         # Check Relitive Error
        if ((iterations > 1) and (midpoint != 0)):
            relError = abs(midpoint - previousMid) / abs(midpoint)

        # Save previous midpoint
        previousMid = midpoint

    # Return Array of values [L, R, midpoint, funcMidpoint, iterations, error]
    return [xLeft, xRight, midpoint, value, iterations, (relError * 100)]
   


def main():
    bisectionList = []

    # Bounding values L and R [1, 2] for -X^2 - 2
    lBound = 1 # From graphical estimation on assignment doc
    rBound = 2  # From graphical estimation on assignment doc

    # Tolerance and max iterations
    tolerance = 1e-5 # 10% relitive error allowed
    maxIterations = 20

    # true: root exists, false: root may exist, none: root is lower bound
    rootSign = checkRootExists(lBound, rBound)


    # If root exists
    if (rootSign < 0):

        # Call bisection method Args (L, R, Tolerance, Max Iterations)
        # Return Array of values [L, R, midpoint, funcMidpoint, iterations, error]
        bisectionList = bisectionMethod(lBound, rBound, tolerance, maxIterations, rootFunc)
        
        # Print Results
        print(f"Bisection Method Results - Bounds: [{lBound}, {rBound}]")    
        print("Root Estimate \t f(Root) \t Iterations \t Relative Error (%)")
        print(round(bisectionList[2], 4), "\t\t", round(bisectionList[3], 4), "\t\t", round(bisectionList[4], 4), "\t\t", round(bisectionList[5], 4))

        # Print final interval and interval size
        finalLeft = bisectionList[0]
        finalRight = bisectionList[1]
        intervalSize = abs(finalRight - finalLeft)
        print(f"Final Interval: [{round(finalLeft, 10)}, {round(finalRight, 10)}]")
        print(f"Final Interval Size: {round(intervalSize, 10)}")

        # Estimate the maximum absolute error when taking the midpoint as the root
        maxAbsError = intervalSize / 2
        print(f"Maximum Absolute Error Estimate: {round(maxAbsError, 10)}")

    else:
        print("Root is: ", lBound)

# Main End

if __name__ == "__main__":
    main()