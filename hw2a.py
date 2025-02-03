#region imports
from math import sqrt, pi, exp
from NumericalMethods import GPDF, Simpson, Probability
#endregion

#region function definitions
def main():
    """
    I want to integrate the Gaussian probability density function between
    a left hand limit = (mean - 5*stDev) to a right hand limit = (c).  Here
    is my step-by-step plan:
    1. Decide mean, stDev, and c and if I want P(x>c) or P(x<c).
    2. Define args tuple and c to be passed to Probability
    3. Pass args, and a callback function (GPDF) to Probability
    4. In probability, pass along GPDF to Simpson along with the appropriate args tuple
    5. Return the required probability from Probability and print to screen.
    :return: Nothing to return, just print results to screen.
    """
    # Test 1: P(x < 105 | N(100, 12.5))
    mean1, stDev1, c1 = 100, 12.5, 105
    p1 = Probability(GPDF, (mean1, stDev1), c1, False)
    print(f"P(x<{c1:.2f}|N({mean1},{stDev1}))={p1:.2f}")

    # Test 2: P(x > μ + 2σ | N(100, 3))
    mean2, stDev2 = 100, 3
    c2 = mean2 + 2 * stDev2
    p2 = Probability(GPDF, (mean2, stDev2), c2, True)
    print(f"P(x>{c2:.2f}|N({mean2},{stDev2}))={p2:.2f}")

#endregion

if __name__ == "__main__":
    main()

