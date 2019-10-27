from scipy.stats import linregress
import math
#import scipy

def stdev(data):
    diffs = 0
    avg = sum(data)/len(data)
    for n in data:
        diffs += (n - avg)**2
    return (diffs/(len(data)-1)) ** 0.5


def math_expectation(data):
    sum_data = sum(data)
    num_data = len(data)
    return sum_data / num_data

def get_coeff_determination(x, y, method):
    if method == 'lin-lin':
        (coefficients, intercept, rvalue, pvalue, stderr) = linregress(x, y)
    elif method == 'lin-log':
        y = math.log2(y)
        (coefficients, intercept, rvalue, pvalue, stderr) = linregress(x, y)
    elif method == 'log-lin':
        x = math.log2(x)
        (coefficients, intercept, rvalue, pvalue, stderr) = linregress(x, y)
    elif method == 'log-log':
        x = math.log2(x)
        y = math.log2(y)
        (coefficients, intercept, rvalue, pvalue, stderr) = linregress(x, y)
    return coefficients, intercept, rvalue, pvalue, stderr