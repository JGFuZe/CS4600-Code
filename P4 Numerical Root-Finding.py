# Defines the function for root finding
# -12 - 21x + 18x^2 - 2.75x^3
def rootFunc(x):
    return ((-12) - (21*x) + (18)*(x**2) - (2.75)*(x**3))


def checkRootExists(lBound, rBound):
    #
    returnVal = 0

    #
    rootLower = (rootFunc(lBound))

    # if f(x) Left Bound Value * Right Bound Value is < 0 Then There is a root
    if ((rootFunc(lBound) * rootFunc(rBound)) < 0):
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
   
def falsePositionMethod(lBound, rBound, tolerance, maxIterations, rootFunc):
    # Initialize Variables
    root = 0
    rootFuncValue = 0
    previousRoot = 0

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

        # Calculate function values at bounds
        funcValRight = rootFunc(xRight)
        funcValLeft = rootFunc(xLeft)

        # Approximate root using false position formula
        root = xRight - ((funcValRight * (xLeft - xRight)) / (funcValLeft - funcValRight))

        # Calculate value at iteration value
        rootFuncValue = rootFunc(root)

        # Update bounds based on sign of function at root
        if ((funcValLeft * rootFuncValue) < 0):
            xRight = root 
        elif ((funcValLeft * rootFuncValue) > 0):
            xLeft = root
        else:
            # Root is found
            break

        # Check Relitive Error
        if ((iterations > 1) and (root != 0)):
            relError = abs(root - previousRoot) / abs(root)

        # Save previous midpoint
        previousRoot = root

    # Return Array of values [L, R, root, rootFunctionValue, iterations, error]
    return [xLeft, xRight, root, rootFuncValue, iterations, (relError * 100)]


def main():
    bisectionList = []
    falsePosList = []

    # Bounding values L and R
    lBound = -1 # From graphical estimation on assignment doc
    rBound = 1  # From graphical estimation on assignment doc

    # Tolerance and max iterations
    tolerance = 0.01 # 10% relitive error allowed
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
        print(round(bisectionList[2], 4), "\t", round(bisectionList[3], 4), "\t", round(bisectionList[4], 4), "\t\t", round(bisectionList[5], 4))


        # Call false position method Args (L, R, Tolerance, Max Iterations)
        # Return Array of values [L, R, midpoint, funcMidpoint, iterations, error]
        falsePosList = falsePositionMethod(lBound, rBound, tolerance, maxIterations, rootFunc)
            # Return Array of values [L, R, midpoint, function(midpoint), error]

        # Print Results
        print(f"\n\nFalse Position Results - Bounds: [{lBound}, {rBound}]")    
        print("Root Estimate \t f(Root) \t Iterations \t Relative Error (%)")
        print(round(falsePosList[2], 4), "\t\t", round(falsePosList[3], 4), "\t", round(falsePosList[4], 4), "\t\t", round(falsePosList[5], 4))

    else: #
        print("Root is: ", lBound)

# Main End

if __name__ == "__main__":
    main()