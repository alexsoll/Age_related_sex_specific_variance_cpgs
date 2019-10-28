import numpy as np
import math
from sklearn.linear_model import LinearRegression

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
        x = np.array(x).reshape(-1, 1)
        y = np.array(y)

        model = LinearRegression()
        model.fit(x, y)
        r_sq = model.score(x, y)
    elif method == 'lin-log':
        y_ = [math.log2(i) for i in y]
        x = np.array(x).reshape(-1, 1)
        y_ = np.array(y_)

        model = LinearRegression()
        model.fit(x, y_)
        r_sq = model.score(x, y)
    elif method == 'log-lin':
        x_ = [math.log2(i) for i in x]
        x_ = np.array(x_).reshape(-1, 1)
        y = np.array(y)

        model = LinearRegression()
        model.fit(x_, y)
        r_sq = model.score(x_, y)
    elif method == 'log-log':
        x_ = [math.log2(i) for i in x]
        y_ = [math.log2(i) for i in y]
        x_ = np.array(x_).reshape(-1, 1)
        y_ = np.array(y_)

        model = LinearRegression()
        model.fit(x_, y_)
        r_sq = model.score(x_, y_)
    return  r_sq