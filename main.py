import argparse
from os import stat_result
#import numpy as np
import processing
import data
#import numpy as np
import stats

def argparser():
    parser = argparse.ArgumentParser(description="Age-Related-sex-specific-variance-cpgs")
    parser.add_argument('--gse', metavar="file", type=str, help="Path to the file containing the GSE database")
    parser.add_argument('-s', '--path-to-save', metavar="file", type=str, help="Path to the file where will the results"
                                                                         "be saved")
    return parser.parse_args()


def test():
    ages = data.get_ages('C:\\Users\\alsol\\Desktop\\scientific_adviser\\attributes_GSE87571.txt')
    with open('C:\\Users\\alsol\\Desktop\\scientific_adviser\\average_beta.txt') as file:
        file.readline()
        '''while True:
            line = file.readline().split()
            if len(line) == 0:
                break
            name = line.pop(0)
            betas = np.asarray([float(beta) for beta in line], dtype=float)
            sigma = stats.stdev(betas)
            avg = stats.math_expectation(betas)
            print(avg + 2 * sigma)
            print(avg - 2 * sigma)'''
        line = file.readline().split()
        name = line.pop(0)
        betas = line
        sigma = stats.stdev(betas)
        avg = stats.math_expectation(betas)
        low, high = processing.get_low_high_betas(betas, ages)
        print(low)
        print(high)
        r2_h = stats.get_coeff_determination(ages, high, 'lin-lin')


def main():
    with open('C:\\Users\\alsol\\Desktop\\scientific_adviser\\average_beta.txt') as file:
        file.readline()
        while True:
            line = file.readline().split()
            if len(line) == 0:
                break
            name = line.pop(0)
            betas = line
            sigma = stats.stdev(betas)
            avg = stats.math_expectation(betas)
            low, high = processing.get_low_high_betas(betas, data.get_ages(
                'C:\\Users\\alsol\\Desktop\\scientific_adviser\\attributes_GSE87571.txt'))




if __name__ == "__main__":
    test()
    #processing.get_low_high_beta()