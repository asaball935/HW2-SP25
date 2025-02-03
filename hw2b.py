#region imports
from NumericalMethods import Secant  # Importing the Secant method from NumericalMethods.py
from math import cos  # Import the cos function from the math module
#endregion

#region function definitions
def fn1(x):
    """
    Function: x - 3 * cos(x) = 0
    We want to find the root of this function.
    """
    return x - 3 * cos(x)


def fn2(x):
    """
    Function: cos(2x) * x^3 = 0
    We want to find the root of this function.
    """
    return cos(2 * x) * x**3


def main():
    """
       fn1:  x - 3cos(x) = 0; with x0=1, x1= 2, maxiter = 5 and xtol = 1e-4
       fn2:  cos(2x) * x^3 = 0; with x0=1, x1= 2, maxiter = 15 and xtol = 1e-8
       fn2:   with x0=1, x1= 2, maxiter = 3 and xtol = 1e-8

       I observe that for functions 2 and 3, the answer should be pi/2 or about 1.57
    :return: nothing, just print results
    """
    # Using Secant method to find roots of fn1 and fn2 with different iterations and tolerances
    r1 = Secant(fn1, 1, 2, 5, 1e-4)  # Call Secant for fn1
    r2 = Secant(fn2, 1, 2, 15, 1e-8)  # Call Secant for fn2 (maxiter=15, xtol=1e-8)
    r3 = Secant(fn2, 1, 2, 3, 1e-8)   # Call Secant for fn2 (maxiter=3, xtol=1e-8)

    # Print the results for each case
    print(f"Root of fn1 = {r1:0.4f}, after {5} iterations")  # since r1 is the final root approximation
    print(f"Root of fn2 (maxiter=15, xtol=1e-8) = {r2:0.4f}, after 15 iterations")
    print(f"Root of fn2 (maxiter=3, xtol=1e-8) = {r3:0.4f}, after 3 iterations")

#endregion

if __name__ == "__main__":
    main()
