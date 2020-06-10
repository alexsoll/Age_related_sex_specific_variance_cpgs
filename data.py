import numpy as np

def get_ages_info(attributes_file, base_name):
    ages = []
    ages_M_indexes = []
    ages_F_indexes = []
    ages_M = []
    ages_F = []
    # ages = np.array([])
    # ages_M_indexes = np.array([])
    # ages_F_indexes = np.array([])
    # ages_M = np.array([])
    # ages_F = np.array([])
    age_indexes_dict_M = {}
    age_indexes_dict_F = {}
    index = 0
    with open(attributes_file, 'r') as ages_file:
        header = ages_file.readline().split()
        age_column_index = [i for i,x in enumerate(header) if 'age' in x][0]
        gender_column_index = [i for i,x in enumerate(header) if 'gender' in x][0]
        while True:
            line = ages_file.readline()
            if line == "":
                break
            line = line.split()
            ages.append(int(line[age_column_index]))
            if line[gender_column_index] == 'M':
                ages_M.append(int(line[age_column_index]))
                if '55763' in base_name and int(line[age_column_index]) < 35:
                    print("Exclude people yanger then 35 years old with index:", index)
                    index += 1
                    continue 
                ages_M_indexes.append(index)
                if line[age_column_index] not in age_indexes_dict_M:
                    age_indexes_dict_M[line[age_column_index]] = [index]
                else:
                    age_indexes_dict_M[line[age_column_index]].append(index)
                index += 1
            if line[gender_column_index] == 'F':
                ages_F.append(int(line[age_column_index]))
                if '55763' in base_name and int(line[age_column_index]) < 35:
                    print("Exclude people yanger then 35 years old with index:", index)
                    index += 1
                    continue    
                ages_F_indexes.append(index)
                if line[age_column_index] not in age_indexes_dict_F:
                    age_indexes_dict_F[line[age_column_index]] = [index]
                else:
                    age_indexes_dict_F[line[age_column_index]].append(index)
                index += 1
    return ages, ages_M, ages_F, ages_M_indexes, ages_F_indexes, age_indexes_dict_M, age_indexes_dict_F



def get_ages(file):
    with open(file, 'r') as ages_file:
        header = ages_file.readline().split()
        age_column_index = [i for i,x in enumerate(header) if 'age' in x][0]
        ages = []
        while True:
            line = ages_file.readline()
            if line == "":
                break
            line = line.split()
            ages.append(int(line[age_column_index]))
    return ages


def get_M_F_ages_indexes(file):
    with open(file, 'r') as ages_file:
        header = ages_file.readline().split()
        age_column_index = [i for i,x in enumerate(header) if 'age' in x][0]
        gender_column_index = [i for i,x in enumerate(header) if 'gender' in x][0]
        ages_M = []
        ages_F = []
        index = 0
        while True:
            line = ages_file.readline()
            if line == "":
                break
            line = line.split()
            if file.find("55763") != -1 and int(line[age_column_index]) < 35:
                    print("Exclude people yanger then 35 years old with index:", index)
                    index += 1
                    continue
            if line[gender_column_index] == 'M':
                ages_M.append(index)
            elif line[gender_column_index] == 'F':
                ages_F.append(index)
            index += 1
    return ages_M, ages_F


def get_age_indexes(file):
    with open(file, 'r') as ages_file:
        header = ages_file.readline().split()
        age_column_index = [i for i,x in enumerate(header) if 'age' in x][0]
        age_indexes = {}
        index = 0
        while True:
            line = ages_file.readline()
            if len(line) == 0:
                break
            info = line.split()
            if file.find("55763") != -1 and int(info[age_column_index]) < 35:
                    print("Exclude people yanger then 35 years old with index:", index)
                    index += 1
                    continue
            if info[2] not in age_indexes:
                age_indexes[info[age_column_index]] = [index]
            else:
                age_indexes[info[age_column_index]].append(index)
            index += 1
    return age_indexes


def get_exclude_cpgs_list(file):
    exclude_list = []
    with open(file) as rfile:
        while True:
            line = rfile.readline()
            if len(line) == 0:
                break
            cpg_name = line.strip()
            exclude_list.append(cpg_name)
    return exclude_list


def get_cpg_in_sex_ch(file):
    exclude_list = []
    with open(file) as rfile:
        while True:
            line = rfile.readline()
            if len(line) == 0:
                break
            cpg_info = line.split("\t")
            if cpg_info[1] == 'X' or cpg_info[1] == 'Y':
                exclude_list.append(cpg_info[0].strip())
    return exclude_list


def get_perc_z(file, perc):
    z = []
    with open(file, 'r') as rfile:
        rfile.readline()
        while True:
            line = rfile.readline()
            if len(line) == 0:
                break
            line = line.strip().split()
            try:
                z.append(float(line[13]))
            except:
                print(line[13])
    z_perc = np.percentile(z, perc)
    return z_perc


def get_best_cpg(file, perc):
    d = 1.5
    best_cpg = []
    with open(file, 'r') as rfile:
        header = rfile.readline()
        header = header.strip().split()
        K_column_index = [i for i,x in enumerate(header) if 'K' == x][0]
        Z_column_index = [i for i,x in enumerate(header) if 'z' == x][0]
        while True:
            line = rfile.readline().strip()
            if len(line) == 0:
                break
            line = line.split("\t")
            if float(line[Z_column_index]) >= float(perc) and float(line[K_column_index]) >= d:
                best_cpg.append(line[0])
    return set(best_cpg)
