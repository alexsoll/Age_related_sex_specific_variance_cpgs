import data
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab
import seaborn as sns
import pandas as pd
import os
from scipy.stats import norm

def get_low_high_betas(betas, ages, age_indexes):
    ages_ = [int(i) for i in list(age_indexes.keys())]
    min_age = min(ages_)
    max_age = max(ages_)
    high_betas = []
    low_betas = []
    vicinity = 4
    using_ages = [i for i in range(min_age, max_age+1, 1)]
    for age in range(min_age, max_age + 1, 1):
        neighboring_betas_ = neighboring_betas(age, ages, betas, age_indexes, vicinity)
        neighboring_betas_ = np.asarray(neighboring_betas_)
        pers_5 = np.percentile(neighboring_betas_, 5)
        pers_95 = np.percentile(neighboring_betas_, 95)
        high_betas.append(pers_95)
        low_betas.append(pers_5)
    return low_betas, high_betas, using_ages

def neighboring_betas(age, ages, betas, age_indexes, vicinity):
    neigh_betas = []
    min_age = min(ages)
    max_age = max(ages)

    l_b = max(min_age, age - vicinity)
    r_b = min(max_age, age + vicinity)
    for i in range(l_b, r_b + 1):
        if str(i) in age_indexes:
            for index in age_indexes[str(i)]:
                neigh_betas.append(betas[index])
    return neigh_betas

def get_age_historam(attributes_file, base_name):
    ages = data.get_ages(attributes_file)
    ages_M_indexes, ages_F_indexes = data.get_M_F_ages_indexes(attributes_file)
    F_ages = []
    M_ages = []
    for i, age in enumerate(ages):
        if i in ages_M_indexes:
            M_ages.append(age)
        else:
            F_ages.append(age)
    plt.title(base_name)
    plt.hist(M_ages, color='deepskyblue', alpha=0.5, bins=35, range=(-5, 105), label='Мужчины')
    plt.hist(F_ages, color='red', alpha=0.5, bins=35, range=(-5, 105), label='Женщины')
    plt.legend()
    plt.ylabel("Количество пациентов", fontsize = 13)
    plt.xlabel('Возраст', fontsize = 13)
    plt.show()


def get_length_sides(betas, ages, age_indexes, best_top_method, best_bottom_method, coef_top, intercept_top,
                     coef_bot, intercept_bot):
    min_age = min(ages)
    max_age = max(ages)
    vicinity = 5
    length = []

    for age in min_age, max_age:
        if best_bottom_method == 'lin-lin':
            low_betas = coef_bot * age + intercept_bot
        elif best_bottom_method == 'lin-log':
            low_betas = coef_bot * age + intercept_bot
            low_betas = math.exp(low_betas)
        elif best_bottom_method == 'log-lin':
            low_betas = coef_bot * math.log(age) + intercept_bot
        elif best_bottom_method == 'log-log':
            low_betas = coef_bot * math.log(age) + intercept_bot
            low_betas = math.exp(low_betas)

        if best_top_method == 'lin-lin':
            high_betas = coef_top * age + intercept_top
        elif best_top_method == 'lin-log':
            high_betas = coef_top * age + intercept_top
            high_betas = math.exp(high_betas)
        elif best_top_method == 'log-lin':
            high_betas = coef_top * math.log(age) + intercept_top
        elif best_top_method == 'log-log':
            high_betas = coef_top * math.log(age) + intercept_top
            high_betas = math.exp(high_betas)

        length.append(high_betas - low_betas)
    return length

def get_z_graphics(bases):
    fig = plt.figure()
    z = []
    colors = ['b', 'r', 'y', 'g']
    perc_75_list = [[], [], []]
    plt.figure(dpi=80)
    for base in list(bases['bases'].keys()):
        with open(bases['bases'][base]['results_table'], 'r') as rfile:
            header = rfile.readline().split()
            Z_column_index = [i for i,x in enumerate(header) if 'z' in x][0]
            while True:
                line = rfile.readline()
                if len(line) == 0:
                    break
                line = line.strip().split()
                z.append(float(line[Z_column_index]))
        perc_75 = np.percentile(z, 75)
        perc_75_list[0].append(base)
        perc_75_list[1].append(perc_75)
        perc_75_list[2].append(colors[0])
        z = sorted(z)
        sns.kdeplot(z, color=colors[0], label=base, alpha=1, linewidth=1.5)
        colors.pop(0)
    plt.ylabel("Плотность распределения")
    plt.xlabel('Коэффициент детерминации')
    plt.vlines(perc_75_list[1], 0, 10, color=["b", 'r', 'y', 'g'])
    plt.legend()
    #plt.show()
    plt.savefig(bases['z_graphs'])
    print("Z graphics was saved into ", bases['z_graphs'])


def get_x_y(ages, method, coefs):
    if method == 'lin-lin':
        betas = [float(coefs[0]) * float(age) + float(coefs[1]) for age in ages]
    if method == 'lin-log':
        betas = [math.exp(float(coefs[0]) * float(age) + float(coefs[1])) for age in ages]
    if method == 'log-lin':
        betas = [(float(coefs[0]) * math.log(float(age)) + float(coefs[1])) for age in ages]
    if method == 'log-log':
        betas = [math.exp(float(coefs[0]) * math.log(float(age)) + float(coefs[1])) for age in ages]
    return ages, betas


def get_graph(bases, base):
    avr_beta_file = bases['bases'][base]['betas']
    attrib_file = bases['bases'][base]['attributes']
    base_name = base
    cpgs_file = bases['cpgs_for_plot']
    annotations = bases['annotations']
    results_file = bases['bases'][base]['results_table']
    cpgs = []
    coefs_list = []
    methods_list = []
    cpg_gene_dict = {}

    cpg_dict = {}

    with open(cpgs_file, 'r') as cpgs_file:
        while True:
            line = cpgs_file.readline()
            if len(line) == 0:
                break
            cpg_name = line.strip()
            cpgs.append(cpg_name)
    cpgs = sorted(cpgs)


    #for cpg in cpgs:
    with open(results_file, 'r') as res_file:
        header = res_file.readline().strip().split()
        K_type_column_index = [i for i,x in enumerate(header) if 'K_type' == x][0]
        while True:
            line = res_file.readline()
            if len(line) == 0:
                break
            line = line.strip().split()
            name = line[0]
            if name in cpgs:
                cpg_dict[name] = {'coefs': [], 'methods': [], 'k_type': '', 'genes': set()}
                cpg_dict[name]['coefs'] = line[5:13]
                cpg_dict[name]['methods'] = line[1:5]
                cpg_dict[name]['k_type'] = line[K_type_column_index]
                coefs = line[5:13]
                methods = line[1:5]
                coefs_list.append(coefs)
                methods_list.append(methods)
    
    #
    #   Find gene for each cpg
    #
    with open(annotations, 'r') as annotation_file:
        annotation_file.readline()
        for line in annotation_file.readlines():
            line = line.strip().split('\t')
            if line[0] in cpgs:
                genes = line[5].strip().split(';')
                cpg_dict[line[0]]['genes'] = set(genes)
    get_cpgs_graphics(avr_beta_file, attrib_file, cpg_dict, base_name, bases)


def get_cpgs_graphics(avr_beta_file, attrib_file, cpgs_dict, base_name, bases):
    ages_M_indexes, ages_F_indexes = data.get_M_F_ages_indexes(attrib_file)
    ages = data.get_ages(attrib_file)
    new_ages = []
    for i, age in enumerate(ages):
        if i in ages_M_indexes or i in ages_F_indexes:
            new_ages.append(age)
    ages = new_ages
    index_cpg = 0
    with open(avr_beta_file, 'r') as rfile:
        rfile.readline()
        graph_index = 1
        while True:
            index_cpg = 0
            line = rfile.readline()
            if len(line) == 0:
                break
            cpg = line.strip().split("\t")
            name = cpg.pop(0)
            if name in list(cpgs_dict.keys()):
                fig = plt.figure()
                for item in list(cpgs_dict.keys()):
                    if item == name:
                        break
                    index_cpg += 1
                cpg = [float(beta) for beta in cpg]
                new_betas = []
                for i, beta_value in enumerate(cpg):
                    if i in ages_M_indexes or i in ages_F_indexes:
                        new_betas.append(beta_value)
                cpg = new_betas

                y_top_bottom_f = []
                y_top_bottom_m = []

                sort_ages = sorted(ages)
                for index, age in enumerate(ages):
                    if index in ages_F_indexes:
                        plt.scatter(age, cpg[index], marker='.', color='r', alpha=0.85, edgecolors='black', linewidth='0.2')
                    else:
                        plt.scatter(age, cpg[index], marker='.', color='deepskyblue', alpha=0.85, edgecolors='black', linewidth='0.2')

                man_y_lim_value = []
                woman_y_lim_value = []
                if cpgs_dict[name]['k_type'] == 'M':
                    index = 4
                    st = []
                    fnsh = []
                    for meth in cpgs_dict[name]['methods'][1], cpgs_dict[name]['methods'][0]:
                        x, y = get_x_y(sort_ages, meth, cpgs_dict[name]['coefs'][index:index+2])
                        y_top_bottom_f.append(y)
                        woman_y_lim_value.append([y[0],y[-1]])
                        fnsh.append(y[-1])
                        st.append(y[0])
                        index += 2
                        plt.plot(x, y, color="red", alpha=0.75, zorder=-1)
                    st = []
                    fnsh = []
                    index = 0
                    for meth in cpgs_dict[name]['methods'][3], cpgs_dict[name]['methods'][2]:
                        x, y = get_x_y(sort_ages, meth, cpgs_dict[name]['coefs'][index:index+2])
                        fnsh.append(y[-1])
                        y_top_bottom_m.append(y)
                        man_y_lim_value.append([y[0],y[-1]])
                        st.append(y[0])
                        index += 2
                        plt.plot(x, y, color="deepskyblue", alpha=0.75,zorder=3)
                else:
                    st = []
                    fnsh = []
                    index = 0
                    for meth in cpgs_dict[name]['methods'][3], cpgs_dict[name]['methods'][2]:
                        x, y = get_x_y(sort_ages, meth, cpgs_dict[name]['coefs'][index:index+2])
                        fnsh.append(y[-1])
                        y_top_bottom_m.append(y)
                        man_y_lim_value.append([y[0],y[-1]])
                        st.append(y[0])
                        index += 2
                        plt.plot(x, y, color="deepskyblue", alpha=0.75, zorder=-1)
                    index = 4
                    st = []
                    fnsh = []
                    for meth in cpgs_dict[name]['methods'][1], cpgs_dict[name]['methods'][0]:
                        x, y = get_x_y(sort_ages, meth, cpgs_dict[name]['coefs'][index:index+2])
                        y_top_bottom_f.append(y)
                        woman_y_lim_value.append([y[0],y[-1]])
                        fnsh.append(y[-1])
                        st.append(y[0])
                        index += 2
                        plt.plot(x, y, color="red", alpha=0.75, zorder=3)
                y_lim_max = max(max(cpg), woman_y_lim_value[0][0], woman_y_lim_value[0][1], man_y_lim_value[0][0], man_y_lim_value[0][1])
                y_lim_min = min(min(cpg), woman_y_lim_value[0][0], woman_y_lim_value[0][1], man_y_lim_value[0][0], man_y_lim_value[0][1])
                title = name
                if list(cpgs_dict[name]['genes'])[0] != '':
                    title += '('
                    for i, item in enumerate(cpgs_dict[name]['genes']):
                        if i != 0:
                            title += ', '
                        title += item
                    title += ')'
                
                if cpgs_dict[name]['k_type'] == 'M':
                    plt.fill_between(x,y_top_bottom_f[0],y_top_bottom_f[1],color='red',alpha=0.5)
                    plt.fill_between(x,y_top_bottom_m[0],y_top_bottom_m[1],color='deepskyblue',alpha=0.65)
                else:
                    plt.fill_between(x,y_top_bottom_m[0],y_top_bottom_m[1],color='deepskyblue',alpha=0.5)
                    plt.fill_between(x,y_top_bottom_f[0],y_top_bottom_f[1],color='red',alpha=0.65)
                plt.xticks(np.arange(0, 130, 50))
                plt.yticks(np.arange(0, 1, 0.1))
                plt.xlim(0, 110)
                plt.ylim(y_lim_min, y_lim_max)
                plt.grid(which='major', axis='x')
                plt.grid(which='major', axis='y')
                plt.title(title)
                plt.xlabel('Возраст')
                plt.ylabel('Уровень метилирования')
                #plt.show()
                fig.savefig(os.path.join(bases['bases'][base_name]['graphs'], name + cpgs_dict[name]['k_type']))
                print('Graphs No. ' + str(graph_index) + ' was saved (' + str(name) + ')')
                graph_index += 1
            else:
                continue