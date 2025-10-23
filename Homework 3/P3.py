# CS4600 - Homework 3 - Problem 3
# Jonah Gallagher


# Bring in NumPy for array math
import numpy as np


# Function to help format matrices nicely
def formatMatrix(matrix):

    # Join each formatted row on its own line
    rowsAsStrings = []
    for row in matrix:
        formattedRow = ", ".join(f"{value: .6e}" for value in row)
        rowsAsStrings.append("[" + formattedRow + "]")
    return "\n".join(rowsAsStrings)

# Function to help format vectors nicely
def formatVector(vector):

    # Wrap vector values in brackets for consistent printing
    formattedValues = ", ".join(f"{value: .6e}" for value in vector)
    return "[" + formattedValues + "]"


# Build the linear system 
def buildReactorSystem():

    # Constants given in the problem statement
    k = 0.1                                                    # min^-1 decay rate
    v1, v2, v3 = 100.0, 50.0, 150.0                            # reactor volumes (L)
    q12, q21, q13, q32, q3out = 5.0, 22.0, 117.0, 7.0, 110.0   # inter-reactor flows (L/min)
    q1_in, q2_in = 100.0, 10.0                                 # inflows (L/min)
    c1_in, c2_in = 10.0, 200.0                                 # feed concentrations (g/L)

    # Assemble A * c = b 
    coefficientMatrix = np.array([
        [-(q12 + q13) - (k * v1), q21, 0.0],  # Reactor 1 row
        [q12, -(q21 + k * v2), q32],          # Reactor 2 row
        [q13, 0.0, -(q32 + q3out + k * v3)]   # Reactor 3 row
    ], dtype=float)

    # Build the right-hand side vector b
    rhsVector = np.array([
        -(q1_in * c1_in),  # Reactor 1 inflow product
        -(q2_in * c2_in),  # Reactor 2 inflow product
        0.0                # Reactor 3 has no external feed
    ], dtype=float)

    # Return tuple of A and b
    return coefficientMatrix, rhsVector
# buildReactorSystem End


# LU factorization and triangular solves ----------------------------------------
def luFactorize(matrix):

    # LU (no pivoting) using compact dot products
    size = matrix.shape[0]  # Get the matrix size
    lower = np.eye(size)    # Initialize L as identity
    upper = np.zeros_like(matrix, dtype=float) # Initialize U as zero matrix

    # Process each pivot column
    for pivot in range(size):

        # Fill the upper matrix entries for the pivot row
        for column in range(pivot, size):
            upper[pivot, column] = matrix[pivot, column] - np.dot(lower[pivot, :pivot], upper[:pivot, column])

        # Fill the lower matrix entries beneath the pivot
        for row in range(pivot + 1, size):

            # Get the current pivot value from upper
            pivotValue = upper[pivot, pivot]

            # Check for zero pivot to avoid division by zero
            if (pivotValue == 0.0):
                raise ZeroDivisionError("Zero pivot encountered during LU factorization.")
            
            # Compute the lower matrix entry
            lower[row, pivot] = (matrix[row, pivot] - np.dot(lower[row, :pivot], upper[:pivot, pivot])) / pivotValue

    # Return L and U
    return lower, upper
# luFactorize End


# Forward substitution function
def forwardSubstitution(lower, rhsVector):

    # Create an empty solution vector for y
    solution = np.zeros_like(rhsVector, dtype=float)

    # March forward through the rows
    for row in range(lower.shape[0]):
        solution[row] = (rhsVector[row] - np.dot(lower[row, :row], solution[:row])) / lower[row, row]

    # Return the intermediate y values
    return solution
# forwardSubstitution End


# Backward substitution function
def backwardSubstitution(upper, rhsVector):

    # Create an empty solution vector for x
    solution = np.zeros_like(rhsVector, dtype=float)

    # March backward from the last row
    for row in range((upper.shape[0] - 1), -1, -1):
        solution[row] = (rhsVector[row] - np.dot(upper[row, (row + 1):], solution[(row + 1):])) / upper[row, row]

    # Return the solved vector
    return solution
# backwardSubstitution End

# Invert matrix using LU factorization
def invertUsingLU(lower, upper):

    # Dimension of the matrix
    size = lower.shape[0]

    # Placeholders for the inverse and identity
    inverseMatrix = np.zeros((size, size), dtype=float)
    identity = np.eye(size)

    # Solve for each column of the inverse
    for column in range(size):
        yValues = forwardSubstitution(lower, identity[:, column])
        inverseMatrix[:, column] = backwardSubstitution(upper, yValues)

    # Return the full inverse matrix
    return inverseMatrix
# invertUsingLU End


# Main Driver -------------------------------------------------------------------
def main():

    # Build the reactor system and factor it
    coefficientMatrix, rhsVector = buildReactorSystem()
    lowerMatrix, upperMatrix = luFactorize(coefficientMatrix)

    # Solve for concentrations and build the inverse
    concentrations = backwardSubstitution(upperMatrix, forwardSubstitution(lowerMatrix, rhsVector))
    inverseMatrix = invertUsingLU(lowerMatrix, upperMatrix)

    # Verify the LU factorization reconstructs A
    reconstructedMatrix = np.matmul(lowerMatrix, upperMatrix)

    # Display all results
    print("Coefficient Matrix [A]:")
    print(formatMatrix(coefficientMatrix))

    print("\nRight-Hand Side {b}:")
    print(formatVector(rhsVector))

    print("\nLower Triangular Matrix [L]:")
    print(formatMatrix(lowerMatrix))

    print("\nUpper Triangular Matrix [U]:")
    print(formatMatrix(upperMatrix))

    print("\nCheck: L @ U reconstruction of [A]:")
    print(formatMatrix(reconstructedMatrix))

    print("\nMatrix Inverse [A^-1]:")
    print(formatMatrix(inverseMatrix))

    print("\nSteady-State Concentrations [c1, c2, c3]^T:")
    print(formatVector(concentrations))

# Main End


if (__name__ == "__main__"):

    main()
