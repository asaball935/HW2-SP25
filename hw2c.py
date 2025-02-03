from NumericalMethods import GaussSeidel
from copy import deepcopy

def main():
    '''
    Solve systems of linear equations using the Gauss-Seidel method
    Example 1: System of equations
    4x + y + z = 7
    x + 3y + 2z = 8
    2x + y + 3z = 9

    Example 2: System of equations
    10x + 2y + 3z = 27
    4x + 15y + 7z = 73
    2x + y + 7z = 29
    '''
    # Example 1: System of equations
    # 3x + y - z = 2
    # x + 4y + z = 12
    # 2x + y + 2z = 10

    # Augmented matrix [A | b]
    Aaug1 = [
        [3, 1, -1, 2],
        [1, 4, 1, 12],
        [2, 1, 2, 10]
    ]

    # Initial guess for the solution vector x
    x1 = [0, 0, 0]

    # Solve using Gauss-Seidel
    solution1 = GaussSeidel(Aaug1, x1, Niter=15)
    print("Solution to system 1:", solution1)

    # Example 2: System of equations
    # 1x - 10y + 2z + 4w = 2
    # 3x + 1y + 4z + 12w = 12
    # 9x + 2y + 3z + 4w = 21
    # -1x + 2y + 7z + 3w = 37

    # Augmented matrix [A | b]
    Aaug2 = [
        [1, -10, 2, 4, 2],
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]
    ]

    # Initial guess for the solution vector x
    x2 = [0, 0, 0, 0]

    # Solve using Gauss-Seidel
    solution2 = GaussSeidel(Aaug2, x2, Niter=15)
    print("Solution to system 2:", solution2)

if __name__ == "__main__":
    main()
