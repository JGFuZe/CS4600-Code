# P6 Newton-Raphson Root-Solving.py

# Defines the function for root finding
# f(x) = x^3 - 6x^2 + 11x - 6.1
def rootFunc(x):
    return ((x**3) - 6*(x**2) + 11*x - 6.1)

# Derivative for Newton-Raphson
def rootDeriv(x):
    return (3*(x**2) - 12*x + 11)

def checkRootExists(lBound, rBound):
    #
    returnVal = 0

    #
    rootLower = (rootFunc(lBound))

    # if f(x) Left Bound Value * Right Bound Value is < 0 Then There is a root
    if ((rootFunc(lBound) * rootFunc(rBound)) < 0):
        returnVal = -1

    elif ((rootLower) == 0):  # Root is the lower bound
        returnVal = 0

    else:  # If not then there is no garenteed root
        returnVal = 1

    return returnVal

def newtonMethod(lBound, rBound, tolerance, maxIterations, rootFunc, rootDeriv):
    # Initialize Variables (match the style of the original methods)
    root = (lBound + rBound) / 2   # use midpoint of the given bounds as the initial guess
    previousRoot = 0
    rootFuncValue = 0

    # iterations counter
    iterations = 0

    # Relative Error (successive-iterate) — initialize similar to other methods
    relError = abs(rBound - lBound)

    # Iterate
    while ((relError > tolerance) and (iterations < maxIterations)):
        # Increment Iterations
        iterations += 1

        f_val = rootFunc(root)
        df_val = rootDeriv(root)

        # Prevent division by zero
        if (df_val == 0):
            break

        # Newton update
        nextRoot = root - (f_val / df_val)

        # Check Relitive Error
        if ((iterations > 0) and (nextRoot != 0)):
            relError = abs(nextRoot - root) / abs(nextRoot)

        # Save previous and update current
        previousRoot = root
        root = nextRoot

    # Final values
    rootFuncValue = rootFunc(root)

    # Return Array of values [L, R, root, rootFunctionValue, iterations, error]
    return [lBound, rBound, root, rootFuncValue, iterations, (relError * 100)]

def main():
    newtonList = []

    # Bounding values L and R (from graphical estimation on assignment doc)
    lBound = 3   # near the largest positive root
    rBound = 4   # near the largest positive root

    # Tolerance and max iterations (keep look-and-feel consistent with your template)
    tolerance = 0.1  # 10% relitive error allowed
    maxIterations = 20

    # true: root exists, false: root may exist, none: root is lower bound
    rootSign = checkRootExists(lBound, rBound)

    # If root exists per the sign check, proceed (keeps the same control flow as your file)
    if (rootSign < 0):
        # Call Newton-Raphson method Args (L, R, Tolerance, Max Iterations)
        # Return Array of values [L, R, root, rootFunctionValue, iterations, error]
        newtonList = newtonMethod(lBound, rBound, tolerance, maxIterations, rootFunc, rootDeriv)

        # Print Results — same format/columns as your original prints
        print(f"Newton-Raphson Results - Bounds: [{lBound}, {rBound}]")
        print("Root Estimate \t f(Root) \t Iterations \t Relative Error (%)")
        print(round(newtonList[2], 4), "\t\t", round(newtonList[3], 4), "\t", round(newtonList[4], 4), "\t\t", round(newtonList[5], 4))

    else:  #
        print("Root is: ", lBound)

# Main End

if __name__ == "__main__":
    main()
