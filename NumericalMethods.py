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
