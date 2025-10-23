# CS4600 - Homework 3 - Problem 2
# Jonah Gallagher

import numpy as np

# Define Pump and Fluid Constants ------------------------------------------------
def pumpFlowRate():
    # Pump discharge flow is 14 liters/min which converts to 0.0002333333 m^3/s
    return (14.0 / 60000.0)
# pumpFlowRate End


def outletPressure():
    # Problem statement sets node 7 (outlet) pressure to 200 kPa
    return 200000.0
# outletPressure End


def formatMatrix(matrix):
    # Build a string where each row appears on a single line for readability
    formattedRows = []
    for row in matrix:
        elements = [f"{value: .6e}" for value in row]
        formattedRows.append("[" + ", ".join(elements) + "]")
    return "\n".join(formattedRows)
# formatMatrix End


def formatVector(vector):
    # Format vector values in a consistent scientific notation
    elements = [f"{value: .6e}" for value in vector]
    return "[" + ", ".join(elements) + "]"
# formatVector End


# Gaussian Elimination with Partial Pivoting -------------------------------------
def gausspivot(coeffMatrix, rhsVector):
    # Convert input arrays to float copies so we can safely modify them
    aMatrix = coeffMatrix.astype(float).copy()
    bVector = rhsVector.astype(float).copy()

    # Grab the matrix size
    size = aMatrix.shape[0]

    # Forward elimination to build upper triangular form
    for pivotIndex in range(size - 1):

        # Locate the row with the largest pivot in the current column
        maxRow = pivotIndex + np.argmax(np.abs(aMatrix[pivotIndex:, pivotIndex]))

        # Swap the rows in both the matrix and RHS vector
        if (maxRow != pivotIndex):
            aMatrix[[pivotIndex, maxRow]] = aMatrix[[maxRow, pivotIndex]]
            bVector[[pivotIndex, maxRow]] = bVector[[maxRow, pivotIndex]]

        # Eliminate the entries below the pivot
        for row in range(pivotIndex + 1, size):
            if (aMatrix[pivotIndex, pivotIndex] == 0):
                continue

            factor = (aMatrix[row, pivotIndex] / aMatrix[pivotIndex, pivotIndex])
            aMatrix[row, pivotIndex:] -= (factor * aMatrix[pivotIndex, pivotIndex:])
            bVector[row] -= (factor * bVector[pivotIndex])

    # Back substitution to solve for the unknowns
    solution = np.zeros(size)
    for row in range((size - 1), -1, -1):
        pivotValue = aMatrix[row, row]
        if (pivotValue == 0):
            solution[row] = 0
            continue

        upperSum = np.dot(aMatrix[row, (row + 1):], solution[(row + 1):])
        solution[row] = ((bVector[row] - upperSum) / pivotValue)

    return solution
# gausspivot End


# Build the 9 x 9 Linear System --------------------------------------------------
def buildFlowSystem():
    # Conductance data (m^3 / (s Pa)) for the pipes
    pipeConductance = {
        "G01": 4.785441e-10,
        "G12": 1.035562e-10,
        "G13": 1.242674e-10,
        "G46": 4.142247e-11,
        "G56": 6.213371e-11,
        "G67": 6.929668e-10
    }

    # Valve coefficients act as conductances (m^3 / (s Pa))
    valveConductance = {
        "C24": 2.00e-09,
        "C35": 2.75e-09
    }

    # Unknown vector order: [P0, P1, P2, P3, P4, P5, P6, Q1, Q2]
    systemMatrix = np.zeros((9, 9))
    rhsVector = np.zeros(9)

    pumpFlow = pumpFlowRate()
    outletP = outletPressure()

    # Helpful indices for readability
    IDX_P0 = 0
    IDX_P1 = 1
    IDX_P2 = 2
    IDX_P3 = 3
    IDX_P4 = 4
    IDX_P5 = 5
    IDX_P6 = 6
    IDX_Q1 = 7
    IDX_Q2 = 8

    # -------- Inlet Pipe --------

    # Equation 1: G01*P0 - G01*P1 = Qpump
    systemMatrix[0, IDX_P0] = pipeConductance["G01"]
    systemMatrix[0, IDX_P1] = (-1 * pipeConductance["G01"])
    rhsVector[0] = pumpFlow


    # -------- Top Branch --------

    # Equation 2: -G12*P1 + G12*P2 + Q1 = 0
    systemMatrix[1, IDX_P1] = (-1 * pipeConductance["G12"])
    systemMatrix[1, IDX_P2] = pipeConductance["G12"]
    systemMatrix[1, IDX_Q1] = 1.0

    # Equation 3: -C24*P2 + C24*P4 + Q1 = 0
    systemMatrix[2, IDX_P2] = (-1 * valveConductance["C24"])
    systemMatrix[2, IDX_P4] = valveConductance["C24"]
    systemMatrix[2, IDX_Q1] = 1.0

    # Equation 4: -G46*P4 + G46*P6 + Q1 = 0
    systemMatrix[3, IDX_P4] = (-1 * pipeConductance["G46"])
    systemMatrix[3, IDX_P6] = pipeConductance["G46"]
    systemMatrix[3, IDX_Q1] = 1.0


    # -------- Bottom Branch --------

    # Equation 5: -G13*P1 + G13*P3 + Q2 = 0
    systemMatrix[4, IDX_P1] = (-1 * pipeConductance["G13"])
    systemMatrix[4, IDX_P3] = pipeConductance["G13"]
    systemMatrix[4, IDX_Q2] = 1.0

    # Equation 6: -C35*P3 + C35*P5 + Q2 = 0
    systemMatrix[5, IDX_P3] = (-1 * valveConductance["C35"])
    systemMatrix[5, IDX_P5] = valveConductance["C35"]
    systemMatrix[5, IDX_Q2] = 1.0

    # Equation 7: -G56*P5 + G56*P6 + Q2 = 0
    systemMatrix[6, IDX_P5] = (-1 * pipeConductance["G56"])
    systemMatrix[6, IDX_P6] = pipeConductance["G56"]
    systemMatrix[6, IDX_Q2] = 1.0


    # -------- After Merge --------
    # Equation 8: G67*P6 = G67*P7 + Qpump
    systemMatrix[7, IDX_P6] = pipeConductance["G67"]
    rhsVector[7] = (pumpFlow + (pipeConductance["G67"] * outletP))

    # -------- Merge continuity --------
    # Equation 9: -Q1 - Q2 = -Qpump
    systemMatrix[8, IDX_Q1] = (-1.0)
    systemMatrix[8, IDX_Q2] = (-1.0)
    rhsVector[8] = (-1 * pumpFlow)

    return (systemMatrix, rhsVector)
# buildFlowSystem End


# Post Processing Utilities ------------------------------------------------------
def litersPerMinute(flowRate):
    # Convert m^3/s to liters/min
    return (flowRate * 60000.0)
# litersPerMinute End


def pascalToPsi(pressure):
    # Convert pressure from Pa to psi using 1 psi = 6895 Pa
    return (pressure / 6895.0)
# pascalToPsi End


def computePipeFlows(solutionVector):
    # Extract pressures and branch flows using the chosen ordering
    p0, p1, p2, p3, p4, p5, p6 = solutionVector[:7]
    q1 = solutionVector[7]
    q2 = solutionVector[8]
    p7 = outletPressure()

    # Build a dictionary of all element flows using the problem orientation
    flows = {
        "G01": pumpFlowRate(),
        "G12": q1,
        "G13": q2,
        "C24": q1,
        "C35": q2,
        "G46": q1,
        "G56": q2,
        "G67": pumpFlowRate()
    }

    # Capture the pressure drops for reporting later
    drops = {
        "G01": (p0 - p1),
        "G12": (p1 - p2),
        "G13": (p1 - p3),
        "C24": (p2 - p4),
        "C35": (p3 - p5),
        "G46": (p4 - p6),
        "G56": (p5 - p6),
        "G67": (p6 - p7)
    }

    return (flows, drops, q1, q2)
# computePipeFlows End


# Main Driver -------------------------------------------------------------------
def main():
    # Build the linear system for the nine unknowns
    coefficientMatrix, rhsVector = buildFlowSystem()

    # Solve using the class gausspivot implementation
    solutionVector = gausspivot(coefficientMatrix, rhsVector)

    # Split the solution for readability
    pNodes = solutionVector[:7]

    # Pre-compute flows and pressure drops for reporting
    elementFlows, elementDrops, q1, q2 = computePipeFlows(solutionVector)

    # Report the matrix equation that was solved
    print("Matrix Equation: [A]{x} = {b}")
    print("Unknown Order: [P0, P1, P2, P3, P4, P5, P6, Q24, Q35]\n")

    print("Coefficient Matrix [A]:")
    print(formatMatrix(coefficientMatrix))
    print("\nRight-Hand Side {b}:")
    print(formatVector(rhsVector))
    print("\nSolution {x}:")
    print(formatVector(solutionVector))
    print("\n---------------------------------------------")

    # Report branch flows in requested units
    print("Branch Flow Summary:")
    print(f"Valve Q24 (Branch 1): {q1:.6e} m^3/s ({litersPerMinute(q1):.4f} L/min)")
    print(f"Valve Q35 (Branch 2): {q2:.6e} m^3/s ({litersPerMinute(q2):.4f} L/min)")
    print(f"Pump Flow Check: {pumpFlowRate():.6e} m^3/s ({litersPerMinute(pumpFlowRate()):.4f} L/min)")
    print("--------------------------------------------------")

    # Report node pressures in both Pascals and psi
    nodeLabels = ["P0", "P1", "P2", "P3", "P4", "P5", "P6"]
    print("Node Pressure Summary:")
    for index, label in enumerate(nodeLabels):
        pressurePa = pNodes[index]
        pressurePsi = pascalToPsi(pressurePa)
        print(f"{label}: {pressurePa:.2f} Pa ({pressurePsi:.2f} psi)")

    outletPa = outletPressure()
    print(f"P7 (Given): {outletPa:.2f} Pa ({pascalToPsi(outletPa):.2f} psi)")
    print("--------------------------------------------------")

    # Detailed element flow report to mirror past assignments
    print("Element Flow Details:")
    for element, flowValue in elementFlows.items():
        dropValue = elementDrops[element]
        print(f"{element}: Flow = {flowValue:.6e} m^3/s ({litersPerMinute(flowValue):.4f} L/min), "
              f"Delta P = {dropValue:.2f} Pa")

# Main End


if (__name__ == "__main__"):
    main()
