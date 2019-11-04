import argparse
from os import stat_result
# import numpy as np
import processing
import data
# import numpy as np
import stats


def argparser():
    parser = argparse.ArgumentParser(description="Age-Related-sex-specific-variance-cpgs")
    parser.add_argument('--gse', metavar="file", type=str, help="Path to the file containing the GSE database")
    parser.add_argument('-s', '--path-to-save', metavar="file", type=str, help="Path to the file where will the results"
                                                                               "be saved")
    return parser.parse_args()


def test():
    ages_M_indexes, ages_F_indexes = data.get_M_F_ages_indexes(
        'C:\\Users\\alsol\\Desktop\\scientific_adviser\\attributes_GSE87571.txt')
    ages = data.get_ages('C:\\Users\\alsol\\Desktop\\scientific_adviser\\attributes_GSE87571.txt')
    with open('C:\\Users\\alsol\\Desktop\\scientific_adviser\\average_beta.txt') as file:
        file.readline()
        line = file.readline().split()
        name = line.pop(0)
        betas = [float(elem) for elem in line]
        low, high, using_ages = processing.get_low_high_betas(betas, ages)
        local_min_r2 = 1
        max_r2 = 0
        for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
            for y in low, high:
                res = stats.get_coeff_determination(using_ages, y, method)
                if local_min_r2 >= res:
                    local_min_r2 = res
            if max_r2 <= local_min_r2:
                max_r2 = local_min_r2

        ########################################################
        local_min_r2_f = 1
        max_r2_f = 0
        ages_F = []
        ages_M = []
        betas_M = []
        betas_F = []
        for i in range(len(ages_F_indexes)):
            ages_F.append(ages[ages_F_indexes[i]])
            betas_F.append(betas[ages_F_indexes[i]])
        for i in range(len(ages_M_indexes)):
            ages_M.append(ages[ages_M_indexes[i]])
            betas_M.append(betas[ages_M_indexes[i]])
        low, high, using_ages = processing.get_low_high_betas(betas_F, ages_F)
        for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
            for y in low, high:
                res_f = stats.get_coeff_determination(using_ages, y, method)
                if local_min_r2_f >= res_f:
                    local_min_r2_f = res_f
            if max_r2_f <= local_min_r2_f:
                max_r2_f = local_min_r2_f
        ####################################################
        local_min_r2_m = 1
        max_r2_m = 0
        low, high, using_ages = processing.get_low_high_betas(betas_M, ages_M)
        for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
            for y in low, high:
                res_m = stats.get_coeff_determination(using_ages, y, method)
                if local_min_r2_m >= res_m:
                    local_min_r2_m = res_m
            if max_r2_m <= local_min_r2_m:
                max_r2_m = local_min_r2_m
        print(max_r2)
        print(max_r2_f)
        print(max_r2_m)
        #######################################################
        length = processing.get_length_sides(betas, ages)
        print(length[0], length[1])


def print_res(table, file):
    dict_keys = ['z', 'r2_M', 'r2_F', 'Start_F', 'Finish_F', 'Start_M', 'Finish_M', 'I_F', 'I_M', 'I_F_type',
                 'I_M_type', 'Y_sex_indep', 'Y_sex_dep']
    print("Write results to a file %s..." % file)
    with open(file, 'w') as wfile:
        line = "CpGs\t"
        for key in dict_keys:
            line += key + '\t'
        line += '\n'
        wfile.write(line)
        line = ''
        for cpg in table.keys():
            line += cpg + '\t'
            for metric in dict_keys:
                line += table[cpg][metric] + '\t'
            line += '\n'
            wfile.write(line)
    print("Write results to a file %s...DONE" % (file))


def main():
    result_table = {}
    count = 1
    ages_M_indexes, ages_F_indexes = data.get_M_F_ages_indexes(
        'C:\\Users\\alsol\\Desktop\\scientific_adviser\\attributes_GSE87571.txt')
    ages = data.get_ages('C:\\Users\\alsol\\Desktop\\scientific_adviser\\attributes_GSE87571.txt')
    with open('C:\\Users\\alsol\\Desktop\\scientific_adviser\\average_beta.txt') as file:
        file.readline()
        while True:
            line = file.readline().split()
            if len(line) == 0:
                break
            name = line.pop(0)
            result_table[name] = {'z': ''}
            betas = [float(elem) for elem in line]
            low, high, using_ages = processing.get_low_high_betas(betas, ages)
            local_min_r2 = 1
            max_r2 = 0
            for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
                for y in low, high:
                    res = stats.get_coeff_determination(using_ages, y, method)
                    if local_min_r2 >= res:
                        local_min_r2 = res
                if max_r2 <= local_min_r2:
                    max_r2 = local_min_r2
            if max_r2 < 0.8:
                continue
            result_table[name]['z'] = max_r2

            ###########################################################

            local_min_r2_f = 1
            max_r2_f = 0
            ages_F = []
            ages_M = []
            betas_M = []
            betas_F = []
            for i in range(len(ages_F_indexes)):
                ages_F.append(ages[ages_F_indexes[i]])
                betas_F.append(betas[ages_F_indexes[i]])
            for i in range(len(ages_M_indexes)):
                ages_M.append(ages[ages_M_indexes[i]])
                betas_M.append(betas[ages_M_indexes[i]])
            low, high, using_ages = processing.get_low_high_betas(betas_F, ages_F)
            for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
                for y in low, high:
                    res_f = stats.get_coeff_determination(using_ages, y, method)
                    if local_min_r2_f >= res_f:
                        local_min_r2_f = res_f
                if max_r2_f <= local_min_r2_f:
                    max_r2_f = local_min_r2_f

            ####################################################

            local_min_r2_m = 1
            max_r2_m = 0
            low, high, using_ages = processing.get_low_high_betas(betas_M, ages_M)
            for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
                for y in low, high:
                    res_m = stats.get_coeff_determination(using_ages, y, method)
                    if local_min_r2_m >= res_m:
                        local_min_r2_m = res_m
                if max_r2_m <= local_min_r2_m:
                    max_r2_m = local_min_r2_m
            #####################################################

            result_table[name]['r2_M'] = max_r2_m
            result_table[name]['r2_F'] = max_r2_f

            ######################################################

            length_M = processing.get_length_sides(betas_M, ages_M)
            length_F = processing.get_length_sides(betas_F, ages_F)

            result_table[name]['Start_F'] = length_F[0]
            result_table[name]['Finish_F'] = length_F[1]
            result_table[name]['Start_M'] = length_M[0]
            result_table[name]['Finish_M'] = length_M[1]

            result_table[name]['I_F'] = max(length_F[0] / length_F[1], length_F[1] / length_F[0])
            result_table[name]['I_M'] = max(length_M[0] / length_M[1], length_M[1] / length_M[0])
            if length_F[0] / length_F[1] > length_F[1] / length_F[0]:
                result_table[name]['I_F_type'] = 0
            else:
                result_table[name]['I_F_type'] = 1

            if length_M[0] / length_M[1] > length_M[1] / length_M[0]:
                result_table[name]['I_M_type'] = 0
            else:
                result_table[name]['I_M_type'] = 1

            #######################################################

            result_table[name]['Y_sex_indep'] = max(result_table[name]['I_M'] / result_table[name]['I_F'],
                                                    result_table[name]['I_F'] / result_table[name]['I_M'])

            #######################################################

            result_table[name]['Y_sex_dep'] = max(max(length_F[0] / length_M[0], length_M[0] / length_F[0]),
                                                  max(length_F[1] / length_M[1], length_M[1] / length_F[1]))
            count += 1
    print_res(result_table, "C:\\Users\\alsol\\Desktop\\scientific_adviser\\Results_table.txt")

if __name__ == "__main__":
    main()
