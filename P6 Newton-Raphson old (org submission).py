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
    rootFuncValue = 0

    # iterations counter
    iterations = 0

    # Relative Error (successive-iterate) — initialize similar to other methods
    relError = abs(rBound - lBound)

    # Iterate
    while ((relError > tolerance) and (iterations < maxIterations)):
        # Increment Iterations
        iterations += 1

        fValue = rootFunc(root)
        dfValue = rootDeriv(root)

        # Prevent division by zero
        if (dfValue == 0):
            break

        # Newton update
        nextRoot = root - (fValue / dfValue)

        # Check Relitive Error
        if ((iterations > 0) and (nextRoot != 0)):
            relError = abs(nextRoot - root) / abs(nextRoot)

        # Save previous and update current
        root = nextRoot

    # Final values
    rootFuncValue = rootFunc(root)

    # Return Array of values [L, R, root, rootFunctionValue, iterations, error]
    return [lBound, rBound, root, rootFuncValue, iterations, (relError * 100)]

def main():
    # List of initial guess intervals
    guessIntervals = [
        [-1.75, 0],
        [0, 1.5], # Worked
        [1.5, 2], # Worked
        [2, 2.5],
        [2.6, 3.5] # Worked
    ]

    tolerance = 0.1  # 10% relative error allowed
    maxIterations = 20

    rootsFound = []

    for interval in guessIntervals:
        lBound, rBound = interval
        rootSign = checkRootExists(lBound, rBound)
        if (rootSign < 0):
            result = newtonMethod(lBound, rBound, tolerance, maxIterations, rootFunc, rootDeriv)
            rootsFound.append(result)
        elif (rootSign == 0):
            # Root is exactly at lBound
            rootsFound.append([lBound, rBound, lBound, rootFunc(lBound), 0, 0])
        # else: no root in this interval

    print("Newton-Raphson Results for All Initial Guesses:")
    print("Interval\tRoot Estimate\tf(Root)\t\tIterations\tRelative Error (%)")
    for res in rootsFound:
        print(f"[{res[0]}, {res[1]}]\t{round(res[2], 4)}\t\t{round(res[3], 4)}\t\t{res[4]}\t\t{round(res[5], 4)}")

# Main End

if __name__ == "__main__":
    main()
