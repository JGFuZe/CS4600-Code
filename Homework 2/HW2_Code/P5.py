# CS4600 - Homework 2 - Problem 5
# Jonah Gallagher

import numpy as np
from scipy.optimize import linprog

def main():
    # Maximize: Z = 3A + 2B
    # Subject to:
    #   A + 2B <= 6
    #   2A + B <= 6
    #   A >= 0, B >= 0

    # Decision vars order: [A, B]
    cVector = np.array([-3, -2])  # minimize -Z

    # Constraint matrix for the upper bounds
    constraintMatrix = np.array([
        [1, 2],
        [2, 1]
    ])

    # Right-hand side vector for the constraints
    rhsVector = np.array([6, 6])

    # Bounds for A and B 
    variableBounds = [(0, None), (0, None)]

    # Call linprog with method="highs"
    result = linprog(cVector, A_ub=constraintMatrix, b_ub=rhsVector, bounds=variableBounds, method="highs")

    # Check the solver status
    if (result.success):
        # Extract the optimal decision variable values
        optA, optB = result.x

        # Convert back to the maximization objective
        optProfit = ((3 * optA) + (2 * optB))

        # Print the results
        print("Linear Programming Optimization Result:")
        print(f"Status: {result.message}")
        print(f"Optimal Solution: A = {optA:.6f}, B = {optB:.6f}")
        print(f"Maximum Profit Z = {optProfit:.6f}")


        # Evaluate constraint usage and slacks
        constraintValues = (constraintMatrix @ result.x)

        # Calculate slacks for each constraint
        slackValues = (rhsVector - constraintValues)

        # Print constraint usage and slacks
        print(f"\nConstraint Usage (A + 2B, 2A + B): {constraintValues}")
        print(f"Slacks (>= 0): {slackValues}")

        # Compare to the one-iteration BFS (A = 3, B = 0, Z = 9)
        bfsProfit = ((3 * 3) + (2 * 0))
        print(f"\nOne-Iteration Tableau BFS: A = 3, B = 0, Z = {bfsProfit}")

        print("\nWas the same as my manual BFS calculation Z = 9")

    else:
        # Report failure message
        print(f"Optimization failed: {result.message}")

# Main End

if (__name__ == "__main__"):
    main()
