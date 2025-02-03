# region imports
import Gauss_Elim as GE  # this is the module from lecture 2 that has useful matrix manipulation functions
from math import sqrt, pi, exp, cos


# endregion

# region function definitions
def Probability(PDF, args, c, GT=True):
    """
    This is the function to calculate the probability that x is >c or <c depending
    on the GT boolean.
    Step 1:  unpack args into mu and stDev
    Step 2:  compute lhl and rhl for Simpson
    Step 3:  package new tuple args1=(mu, stDev, lhl, rhl) to be passed to Simpson
    Step 4:  call Simpson with GNPDF and args1
    Step 5:  return probability
    :param PDF: the probability density function to be integrated
    :param args: a tuple with (mean, standard deviation)
    :param c: value for which we ask the probability question
    :param GT: boolean deciding if we want probability x>c (True) or x<c (False)
    :return: probability value
    """
    mu, sig = args

    # Define the lower and upper bounds for integration
    if GT:
        lhl = c  # left-hand limit is c for x > c
        rhl = mu + 5 * sig  # right-hand limit is μ + 5σ (effectively infinity)
    else:
        lhl = mu - 5 * sig  # left-hand limit is μ - 5σ
        rhl = c  # right-hand limit is c for x < c

    # Perform Simpson's 1/3 rule integration
    p = Simpson(PDF, (mu, sig, lhl, rhl))

    return p


def GPDF(args):
    """
    Here is where I will define the Gaussian probability density function.
    This requires knowing the population mean and standard deviation.
    To compute the GPDF at any value of x, I just need to compute as stated
    in the homework assignment.
    Step 1:  unpack the args tuple into variables called: x, mu, stDev
    Step 2:  compute GPDF value at x
    Step 3:  return value
    :param args: (x, mean, standard deviation)  tuple in that order
    :return: value of GPDF at the desired x
    """
    # Step 1: unpack args
    x, mu, sig = args
    # Step 2: compute GPDF at x
    fx = (1 / (sig * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sig) ** 2)
    # Step 3: return value
    return fx


def Simpson(fn, args, N=100):
    """
    This executes the Simpson 1/3 rule for numerical integration.

    :param fn: The function to integrate (PDF in this case).
    :param args: A tuple containing (mean, stDev, lhl, rhl).
    :param N: The number of subintervals to divide the integration range into.
    :return: The estimated area under the function (i.e., the probability).
    """
    mu, sig, lhl, rhl = args

    # Ensure an even number of intervals
    if N % 2 == 1:
        N += 1

    # Calculate the step size
    h = (rhl - lhl) / N

    # Simpson's 1/3 rule requires evaluating at evenly spaced points
    x_vals = [lhl + i * h for i in range(N + 1)]
    f_vals = [fn((x, mu, sig)) for x in x_vals]  # Evaluate the PDF at each point

    # Apply Simpson's 1/3 rule
    integral = f_vals[0] + f_vals[-1]  # First and last terms
    for i in range(1, N):
        if i % 2 == 0:
            integral += 2 * f_vals[i]  # Even terms multiplied by 2
        else:
            integral += 4 * f_vals[i]  # Odd terms multiplied by 4

    # Multiply by h/3 for Simpson's rule
    integral *= h / 3

    return integral


def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    Implements the Secant method to find the root of a function.

    :param fcn: The function for which we want to find the root.
    :param x0: First guess of the root.
    :param x1: Second guess of the root.
    :param maxiter: Maximum number of iterations before stopping.
    :param xtol: Tolerance level for stopping criterion based on the difference between successive approximations.

    :return: The final estimate of the root (most recent new x value).
    """

    # Iterate up to maxiter times
    for i in range(maxiter):
        # Evaluate function at the two guesses
        f0 = fcn(x0)
        f1 = fcn(x1)

        # Check if we are close enough to the root
        if abs(f1 - f0) < 1e-12:  # Avoid division by zero in case the difference between f1 and f0 is too small
            print("Warning: Small function difference. Possible division by zero.")
            return None

        # Compute the next approximation using the Secant formula
        x_new = x1 - f1 * (x1 - x0) / (f1 - f0)

        # Check the stopping criteria: if the change in x is smaller than xtol, stop the iteration
        if abs(x_new - x1) < xtol:
            print(f"Converged in {i + 1} iterations.")
            return x_new

        # Update x0 and x1 for the next iteration
        x0, x1 = x1, x_new

    # If the loop ends without convergence, print a message and return the last approximation
    print(f"Max iterations reached. Returning last estimate after {maxiter} iterations.")
    return x1


def GaussSeidel(Aaug, x, Niter=15):
    """
    Solve the system of linear equations Ax = b using the Gauss-Seidel method.

    :param Aaug: Augmented matrix [A | b] with N rows and N+1 columns (N equations, N+1 variables).
    :param x: Initial guess for the solution (vector of size N).
    :param Niter: Number of iterations to compute (default is 15).
    :return: The final solution vector after Niter iterations.
    """
    # Extract matrix A and vector b from augmented matrix Aaug
    N = len(Aaug)  # Number of equations (rows)
    A = [row[:-1] for row in Aaug]  # Extract the coefficient matrix A (remove last column)
    b = [row[-1] for row in Aaug]  # Extract the right-hand side vector b (last column)

    # Iterate for the specified number of iterations
    for k in range(Niter):
        # Make a copy of the current solution vector to update in the same iteration
        x_new = x.copy()

        for i in range(N):
            # Sum over all elements except the diagonal element (i-th row)
            sigma = sum(A[i][j] * x_new[j] for j in range(N) if j != i)

            # Update the current solution using the Gauss-Seidel formula
            x_new[i] = (b[i] - sigma) / A[i][i]

        # Update the solution vector x with the new estimates
        x = x_new

    return x


def main():
    '''
    This is a function I created for testing the numerical methods locally.
    :return: None
    '''
    # region testing GPDF
    fx = GPDF((0, 0, 1))
    print("{:0.5f}".format(fx))  # Does this match the expected value?
    # endregion

    # region testing Simpson
    p = Simpson(GPDF, (0, 1, -5, 0))  # should return 0.5
    print("p={:0.5f}".format(p))  # Does this match the expected value?
    # endregion

    # region testing Probability
    p1 = Probability(GPDF, (0, 1), 0, True)
    print("p1={:0.5f}".format(p1))  # Should be approximately 0.5 for P(x > 0) with standard normal distribution
    p2 = Probability(GPDF, (0, 1), 0, False)
    print("p2={:0.5f}".format(p2))  # Should be approximately 0.5 for P(x < 0) with standard normal distribution
    # endregion


# endregion

# region function calls
if __name__ == '__main__':
    main()
# endregion