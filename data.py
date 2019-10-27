#import numpy as np


def get_ages(file):
    with open(file, 'r') as ages_file:
        ages_file.readline()
        ages = []
        while True:
            line = ages_file.readline()
            if line == "":
                break
            ages.append(int(line.split()[3]))
    return ages

def get_M_F_ages_indexes(file):
    with open(file, 'r') as ages_file:
        ages_file.readline()
        ages_M = []
        ages_F = []
        index = 0
        while True:
            line = ages_file.readline()
            if line == "":
                break
            line = line.split(" ")
            if line[2] == 'M':
                ages_M.append(index)
            elif line[2] == 'F':
                ages_F.append(index)
            index += 1
    return ages_M, ages_F


def neighboring_betas(age, ages, betas, vicinity):
    neigh_betas = []
    for index, item in enumerate(ages):
        if item >= age - vicinity and item <= age + vicinity:
            neigh_betas.append(betas[index])
    return neigh_betas