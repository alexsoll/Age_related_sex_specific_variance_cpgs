import argparse
import numpy as np
import processing
import data
import os
# import numpy as np
import stats
import threading
import time
from os import stat_result


def argparser():
    parser = argparse.ArgumentParser(description="Age-Related-sex-specific-variance-cpgs")
    parser.add_argument('-r', '--gse-root', metavar="dir", type=str, help="Path to the root dir", required=True)
    parser.add_argument('-g', '--gse', metavar="dir", type=str, help="The names of GSE database \
                        dirs in root directory (example: gse40279,gse57871,EPIC)", required=True)
    parser.add_argument('-a', '--action', type=str, help="Use: \'calculation\', \'graphs\' or \'all\'", required=True)
    return parser.parse_args()


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
        for cpg in table:
            line += cpg + '\t'
            for metric in dict_keys:
                line += str(table[cpg][metric]) + '\t'
            line += '\n'
            wfile.write(line)
    print("Write results to a file %s...DONE" % (file))


def calculation(bases, base_name):
    exclude_cpgs = data.get_exclude_cpgs_list(bases['bases'][base_name]['bad_cpgs'])
    cpg_in_sex_ch = data.get_cpg_in_sex_ch(bases['annotations'])

    result_table = {}

    #
    #   Get all ages info
    #
    ages, ages_M, ages_F, ages_M_indexes, ages_F_indexes, age_indexes_M, age_indexes_F = data.get_ages_info(bases['bases'][base_name]['attributes'], base_name) 

    dict_keys = ['best_top_method_f', 'best_bottom_method_f', 'best_top_method_m', 'best_bottom_method_m',
                 'coef_bot_m', 'intercept_bot_m', 'coef_top_m', 'intercept_top_m',
                 'coef_bot_f', 'intercept_bot_f', 'coef_top_f', 'intercept_top_f',
                 'z', 'Start_F', 'Finish_F', 'Start_M', 'Finish_M',
                 'I_F', 'I_M', 'I_F_type', 'I_M_type', 'Y', 'K', 'K_type']

    #
    #   Printing a header to a file with results
    #
    with open(bases['bases'][base_name]['results_table'], 'a') as wfile:
        wline = "CpGs\t"
        for key in dict_keys:
            wline += key + '\t'
        wline += '\n'
        wfile.write(wline)

    #
    #   Main cycle with calculations
    #   
    with open(bases['bases'][base_name]['betas'], 'r') as file:
        file.readline()
        while True:
            wline = ''
            line = file.readline().split()
            if len(line) == 0:
                break
            name = line.pop(0)
            if name in exclude_cpgs or name in cpg_in_sex_ch:
                continue
            betas = [float(elem) for elem in line]
            ###########################################################

            result_table[name] = {'z': ''}
            result_table[name]['best_top_method_f'] = ''
            result_table[name]['best_bottom_method_f'] = ''
            result_table[name]['best_top_method_m'] = ''
            result_table[name]['best_bottom_method_m'] = ''
            result_table[name]['coef_bot_m'] = ''
            result_table[name]['intercept_bot_m'] = ''
            result_table[name]['coef_bot_f'] = ''
            result_table[name]['intercept_bot_f'] = ''
            result_table[name]['coef_top_m'] = ''
            result_table[name]['intercept_top_m'] = ''
            result_table[name]['coef_top_f'] = ''
            result_table[name]['intercept_top_f'] = ''

            local_max_r2_f_bottom = -np.inf
            local_max_r2_f_top = -np.inf

            betas_M = []
            betas_F = []

            for index in ages_F_indexes:
                betas_F.append(betas[index])
            for index in ages_M_indexes:
                betas_M.append(betas[index])

            ###########################################################

            low, high, using_ages = processing.get_low_high_betas(betas, ages_F, age_indexes_F)
            for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
                res_f = stats.get_coeff_determination(using_ages, low, method)
                if type(res_f) != list:
                    continue
                else:
                    if local_max_r2_f_bottom <= res_f[0]:
                        local_max_r2_f_bottom = res_f[0]
                        result_table[name]['best_bottom_method_f'] = method
                        result_table[name]['coef_bot_f'] = res_f[1][0]
                        result_table[name]['intercept_bot_f'] = res_f[2]

            for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
                res_f = stats.get_coeff_determination(using_ages, high, method)
                if type(res_f) != list:
                    continue
                else:
                    if local_max_r2_f_top <= res_f[0]:
                        local_max_r2_f_top = res_f[0]
                        result_table[name]['best_top_method_f'] = method
                        result_table[name]['coef_top_f'] = res_f[1][0]
                        result_table[name]['intercept_top_f'] = res_f[2]

            r2_best_f = min(local_max_r2_f_bottom, local_max_r2_f_top)

            ####################################################

            local_max_r2_m_bottom = -np.inf
            local_max_r2_m_top = -np.inf

            low, high, using_ages = processing.get_low_high_betas(betas, ages_M, age_indexes_M)

            for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
                res_m = stats.get_coeff_determination(using_ages, low, method)
                if type(res_m) != list:
                    continue
                else:
                    if local_max_r2_m_bottom <= res_m[0]:
                        local_max_r2_m_bottom = res_m[0]
                        result_table[name]['best_bottom_method_m'] = method
                        result_table[name]['coef_bot_m'] = res_m[1][0]
                        result_table[name]['intercept_bot_m'] = res_m[2]

            for method in 'lin-lin', 'lin-log', 'log-lin', 'log-log':
                res_m = stats.get_coeff_determination(using_ages, high, method)
                if type(res_m) != list:
                    continue
                else:
                    if local_max_r2_m_top <= res_m[0]:
                        local_max_r2_m_top = res_m[0]
                        result_table[name]['best_top_method_m'] = method
                        result_table[name]['coef_top_m'] = res_m[1][0]
                        result_table[name]['intercept_top_m'] = res_m[2]

            r2_best_m = min(local_max_r2_f_bottom, local_max_r2_f_top)

            #####################################################

            result_table[name]['z'] = min(r2_best_m, r2_best_f)

            wline += name + '\t' + str(result_table[name]['best_top_method_f']) \
                     + '\t' + str(result_table[name]['best_bottom_method_f']) \
                     + '\t' + str(result_table[name]['best_top_method_m']) \
                     + '\t' + str(result_table[name]['best_bottom_method_m']) \
                     + '\t' + str(result_table[name]['coef_bot_m']) \
                     + '\t' + str(result_table[name]['intercept_bot_m']) \
                     + '\t' + str(result_table[name]['coef_top_m']) \
                     + '\t' + str(result_table[name]['intercept_top_m']) \
                     + '\t' + str(result_table[name]['coef_bot_f']) \
                     + '\t' + str(result_table[name]['intercept_bot_f']) \
                     + '\t' + str(result_table[name]['coef_top_f']) \
                     + '\t' + str(result_table[name]['intercept_top_f']) \
                     + '\t' + str(result_table[name]['z']) + '\t'

            ######################################################

            length_M = processing.get_length_sides(betas, ages_M, age_indexes_M,
                                                   result_table[name]['best_top_method_m'],
                                                   result_table[name]['best_bottom_method_m'],
                                                   result_table[name]['coef_top_m'],
                                                   result_table[name]['intercept_top_m'],
                                                   result_table[name]['coef_bot_m'],
                                                   result_table[name]['intercept_bot_m'])
            length_F = processing.get_length_sides(betas, ages_F, age_indexes_F,
                                                   result_table[name]['best_top_method_f'],
                                                   result_table[name]['best_bottom_method_f'],
                                                   result_table[name]['coef_top_f'],
                                                   result_table[name]['intercept_top_f'],
                                                   result_table[name]['coef_bot_f'],
                                                   result_table[name]['intercept_bot_f'])

            result_table[name]['Start_F'] = length_F[0]
            result_table[name]['Finish_F'] = length_F[1]
            result_table[name]['Start_M'] = length_M[0]
            result_table[name]['Finish_M'] = length_M[1]

            wline += str(result_table[name]['Start_F']) + '\t' + str(result_table[name]['Finish_F']) \
                     + '\t' + str(result_table[name]['Start_M']) + '\t' + str(result_table[name]['Finish_M']) + '\t'

            ######################################################

            result_table[name]['I_F'] = max(result_table[name]['Start_F'] / result_table[name]['Finish_F'],
                                            result_table[name]['Finish_F'] / result_table[name]['Start_F'])

            result_table[name]['I_M'] = max(result_table[name]['Start_M'] / result_table[name]['Finish_M'],
                                            result_table[name]['Finish_M'] / result_table[name]['Start_M'])

            wline += str(result_table[name]['I_F']) + '\t' + str(result_table[name]['I_M']) + '\t'

            ######################################################

            if length_F[0] / length_F[1] > length_F[1] / length_F[0]:
                result_table[name]['I_F_type'] = 0
            else:
                result_table[name]['I_F_type'] = 1

            if length_M[0] / length_M[1] > length_M[1] / length_M[0]:
                result_table[name]['I_M_type'] = 0
            else:
                result_table[name]['I_M_type'] = 1

            wline += str(result_table[name]['I_F_type']) + '\t' + str(result_table[name]['I_M_type']) + '\t'

            #######################################################

            result_table[name]['Y'] = max(result_table[name]['I_F'], result_table[name]['I_M'])

            wline += str(result_table[name]['Y']) + '\t'

            #######################################################

            result_table[name]['K'] = max(result_table[name]['I_F'] / result_table[name]['I_M'], 
                                            result_table[name]['I_M'] / result_table[name]['I_F'])

            if result_table[name]['I_F'] / result_table[name]['I_M'] > result_table[name]['I_M'] / result_table[name]['I_F']:
                result_table[name]['K_type'] = 'F'
            else:
                result_table[name]['K_type'] = 'M'

            wline += str(result_table[name]['K']) + '\t' + str(result_table[name]['K_type']) + '\n'

            #######################################################

            with open(bases['bases'][base_name]['results_table'], 'a') as wfile:
                wfile.write(wline)


def bases_intersection(bases):
    intersection_total = set()
    for base in list(bases['bases'].keys()):
        base_percentile = data.get_perc_z(bases['bases'][base]['results_table'], 75)
        bases['bases'][base]['best_cpg'] = data.get_best_cpg(bases['bases'][base]['results_table'], base_percentile)
        print("INFO: count of best cpgs for ", base, " : ", len(bases['bases'][base]['best_cpg']))
        print("INFO: percentile ", base, " : ", base_percentile)
        if len(intersection_total) == 0:
            intersection_total = bases['bases'][base]['best_cpg']
        else:
            intersection_total = set.intersection(intersection_total, bases['bases'][base]['best_cpg'])
    print("INFO: number of cpgs in intersection_total = " + str(len(intersection_total)))

    with open(bases['intersection'], 'w') as wfile:
        for item in intersection_total:
            wfile.write(str(item) + '\n')


if __name__ == "__main__":
    bases = {}

    args = argparser()
    bases['root_dir'] = args.gse_root
    bases['bases'] = {}

    gse_names = args.gse.split(',')
    for name in gse_names:
        attributes_file = [os.path.join(bases['root_dir'], name, f) for f in os.listdir(os.path.join(bases['root_dir'], name)) if f.find('attributes') != -1][0]
        betas_file = [os.path.join(bases['root_dir'], name, f) for f in os.listdir(os.path.join(bases['root_dir'], name)) if f.find('betas') != -1][0]
        bad_cpgs_file = [os.path.join(bases['root_dir'], name, f) for f in os.listdir(os.path.join(bases['root_dir'], name)) if f.find('bad_cpg') != -1][0]
        bases['bases'][name] = {
                                'attributes': attributes_file, 'betas': betas_file, 'bad_cpgs': bad_cpgs_file, 
                                'results_table': os.path.join(bases['root_dir'], name, 'Result_table_' + name + '.txt'),
                                'best_cpg' : [], 'graphs': os.path.join(bases['root_dir'], name, 'graphs')
                                }
    bases['annotations'] = [os.path.join(bases['root_dir'], f) for f in os.listdir(bases['root_dir']) if f.find('annotations') != -1][0]
    bases['z_graphs'] = os.path.join(bases['root_dir'], 'Common', 'total_Z_graphs')
    bases['intersection'] = os.path.join(bases['root_dir'], 'Common', 'total_intersection.txt')
    bases['cpgs_for_plot'] = os.path.join(bases['root_dir'], 'Common', 'cpgs_for_plot.txt')

    if args.action == 'calculate':
        for base_name in bases['bases']:
            calculation(bases, base_name)
            print('The base ', base_name, ' was calculated')
    for base in list(bases['bases'].keys()):
        if not os.path.exists(bases['bases'][base]['results_table']):
            print("ERROR: ", bases['bases'][base]['results_table'], " doesn't exist!")
            exit(1)
    if args.action == 'get-best-cpgs':
        bases_intersection(bases)
    if args.action == 'graphs':
        processing.get_z_graphics(bases)
        bases['cpgs_for_plot'] = bases['intersection']
        for base in list(bases['bases'].keys()):
            processing.get_graph(bases, base)       
    #processing.get_graph(base_GSE55763)
    #processing.get_graph(base_GSE55763, base_GSE87571, base_GSE40279, base_EPIC)
    #processing.get_all_gene_list()
    #processing.get_age_historam(bases_root + 'attributes_gse55763.txt', 'GSE55763')