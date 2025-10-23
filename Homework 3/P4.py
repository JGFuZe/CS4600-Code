# CS4600 - Homework 3 - Problem 4
# Jonah Gallagher


# Import NumPy for array handling
import numpy as np


# Helper to format matrices on separate lines -----------------------------------
def formatMatrix(matrix):

    # Store each row string here
    rowStrings = []

    # Build every row string with scientific notation values
    for row in matrix:
        rowStrings.append("[" + ", ".join(f"{value: .6e}" for value in row) + "]")

    # Join the rows with newline characters
    return "\n".join(rowStrings)


# Helper to format vectors on a single line -------------------------------------
def formatVector(vector):

    # Build the contents inside the brackets
    innerText = ", ".join(f"{value: .6e}" for value in vector)

    # Wrap with brackets and return the string
    return "[" + innerText + "]"


# Build the ventilation system matrices -----------------------------------------
def buildRoomSystem(room4Load):

    # Define the airflow connections (m^3/hr)
    q13 = 50.0
    q31 = 50.0
    q23 = 50.0
    q32 = 50.0
    q34 = 90.0
    q43 = 90.0
    q3_out = 50.0

    # Define the external supply rates (m^3/hr)
    q1_in = 150.0
    q2_in = 50.0

    # Define the external supply concentrations (mg/m^3)
    c1_in = 1.0
    c2_in = 40.0

    # Build the coefficient matrix from the mass balances
    coefficientMatrix = np.array([
        [-(q13), 0.0, q31, 0.0],
        [0.0, -(q23), q32, 0.0],
        [q13, q23, -(q31 + q32 + q34 + q3_out), q43],
        [0.0, 0.0, q34, -(q43)]
    ], dtype=float)

    # Build the right-hand side vector using the supplies and load
    rhsVector = np.array([
        -(q1_in * c1_in),
        -(q2_in * c2_in),
        0.0,
        -room4Load
    ], dtype=float)

    # Return the matrix and vector
    return coefficientMatrix, rhsVector
# buildRoomSystem End


# Perform LU factorization (Doolittle) -------------------------------------------
def luFactorize(matrix):

    # Determine the matrix size
    size = matrix.shape[0]

    # Initialize L as identity and U as zeros
    lower = np.eye(size)
    upper = np.zeros_like(matrix, dtype=float)

    # Iterate across each pivot column
    for pivot in range(size):

        # Fill row of U
        for column in range(pivot, size):
            upper[pivot, column] = matrix[pivot, column] - np.dot(lower[pivot, :pivot], upper[:pivot, column])

        # Fill column of L
        for row in range(pivot + 1, size):
            pivotValue = upper[pivot, pivot]
            if (pivotValue == 0.0):
                raise ZeroDivisionError("Zero pivot encountered during LU factorization.")
            lower[row, pivot] = (matrix[row, pivot] - np.dot(lower[row, :pivot], upper[:pivot, pivot])) / pivotValue

    # Return the factorization
    return lower, upper
# luFactorize End


# Forward substitution for L * y = b ---------------------------------------------
def forwardSubstitution(lower, rhsVector):

    # Initialize the solution vector
    solution = np.zeros_like(rhsVector, dtype=float)

    # Sweep forward across rows
    for row in range(lower.shape[0]):
        solution[row] = (rhsVector[row] - np.dot(lower[row, :row], solution[:row])) / lower[row, row]

    # Return the intermediate vector y
    return solution
# forwardSubstitution End


# Backward substitution for U * x = y -------------------------------------------
def backwardSubstitution(upper, rhsVector):

    # Initialize the solution vector
    solution = np.zeros_like(rhsVector, dtype=float)

    # Sweep backward from the bottom row
    for row in range((upper.shape[0] - 1), -1, -1):
        solution[row] = (rhsVector[row] - np.dot(upper[row, (row + 1):], solution[(row + 1):])) / upper[row, row]

    # Return the final vector x
    return solution
# backwardSubstitution End


# Invert a matrix using its LU factors ------------------------------------------
def invertUsingLU(lower, upper):

    # Get the matrix size
    size = lower.shape[0]

    # Prepare the inverse matrix container
    inverseMatrix = np.zeros((size, size), dtype=float)

    # Identity matrix supplies e_i columns
    identity = np.eye(size)

    # Solve for each column of the inverse
    for column in range(size):
        yValues = forwardSubstitution(lower, identity[:, column])
        inverseMatrix[:, column] = backwardSubstitution(upper, yValues)

    # Return the computed inverse
    return inverseMatrix
# invertUsingLU End


# Main Driver -------------------------------------------------------------------
def main():

    # Set the baseline room 4 load (mg/hr)
    baselineLoad = 5000.0

    # Build the system for the baseline load
    coefficientMatrix, rhsVector = buildRoomSystem(baselineLoad)

    # Compute the LU factorization of the coefficient matrix
    lowerMatrix, upperMatrix = luFactorize(coefficientMatrix)

    # Solve the system for current concentrations
    concentrations = backwardSubstitution(upperMatrix, forwardSubstitution(lowerMatrix, rhsVector))

    # Compute the inverse matrix using LU factors
    inverseMatrix = invertUsingLU(lowerMatrix, upperMatrix)

    # Reconstruct the coefficient matrix to verify LU accuracy
    reconstructedMatrix = np.matmul(lowerMatrix, upperMatrix)

    # Compute the constant contribution to room 2 concentration
    constantContribution = np.dot(inverseMatrix[1, :3], rhsVector[:3])

    # Extract the sensitivity of room 2 concentration to the room 4 load
    loadCoefficient = (-1.0) * inverseMatrix[1, 3]

    # Target concentration for room 2 (mg/m^3)
    targetRoom2 = 20.0

    # Solve for the room 4 load required to hit the target
    requiredLoad = (targetRoom2 - constantContribution) / loadCoefficient

    # Compute the reduction needed from the baseline load
    loadReduction = baselineLoad - requiredLoad

    # Report the coefficient matrix
    print("Coefficient Matrix [A]:")
    print(formatMatrix(coefficientMatrix))

    # Report the right-hand side vector
    print("\nRight-Hand Side {b}:")
    print(formatVector(rhsVector))

    # Report the LU factors and reconstruction
    print("\nLower Triangular Matrix [L]:")
    print(formatMatrix(lowerMatrix))
    print("\nUpper Triangular Matrix [U]:")
    print(formatMatrix(upperMatrix))
    print("\nCheck: L @ U reconstruction of [A]:")
    print(formatMatrix(reconstructedMatrix))

    # Report the inverse matrix
    print("\nMatrix Inverse [A^-1]:")
    print(formatMatrix(inverseMatrix))

    # Report the baseline steady-state concentrations
    print("\nSteady-State Concentrations [c1, c2, c3, c4]^T:")
    print(formatVector(concentrations))

    # Report the required load and reduction to reach the target
    print(f"\nLoad needed for c2 = {targetRoom2} mg/m^3: {requiredLoad:.2f} mg/hr")
    print(f"Required load reduction from baseline: {loadReduction:.2f} mg/hr")

# Main End


if (__name__ == "__main__"):

    main()
