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
    #model = LinearRegression()
    if method == 'lin-lin':
        x = np.array(x).reshape(-1, 1)
        y = np.array(y)

        model = LinearRegression()
        model.fit(x, y)
        r_sq = model.score(x, y)
    elif method == 'lin-log':
        try:
            y_ = [math.log(i) for i in y]
            x = np.array(x).reshape(-1, 1)
            y_ = np.array(y_)

            model = LinearRegression()
            model.fit(x, y_)
            r_sq = model.score(x, y)
        except:
            r_sq = 'NA'
    elif method == 'log-lin':
        try:
            x_ = [math.log(i) for i in x]
            x_ = np.array(x_).reshape(-1, 1)
            y = np.array(y)

            model = LinearRegression()
            model.fit(x_, y)
            r_sq = model.score(x_, y)
        except:
            r_sq = 'NA'
    elif method == 'log-log':
        try:
            x_ = [math.log(i) for i in x]
            y_ = [math.log(i) for i in y]
            x_ = np.array(x_).reshape(-1, 1)
            y_ = np.array(y_)

            model = LinearRegression()
            model.fit(x_, y_)
            r_sq = model.score(x_, y_)
        except:
            r_sq = 'NA'
    if r_sq == 'NA':
        return r_sq
    else:
        return [r_sq, model.coef_, model.intercept_]