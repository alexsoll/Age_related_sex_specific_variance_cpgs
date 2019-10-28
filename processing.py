import stats
import data

def get_low_high_betas(betas, ages):
    min_age = min(ages)
    max_age = max(ages)
    high_betas = []
    low_betas = []
    vicinity = 5
    print("INFO: Getting upper and lower beta values...")
    using_ages = [i for i in range(min_age, max_age+1, 1)]
    for age in range(min_age, max_age+1, 1):
        neighboring_betas = data.neighboring_betas(age, ages, betas, vicinity)
        sigma = stats.stdev(neighboring_betas)
        avg = stats.math_expectation(neighboring_betas)
        high_betas.append(avg + 2 * sigma)
        low_betas.append(avg - 2 * sigma)
    print("INFO: Getting upper and lower beta values...DONE")
    return low_betas, high_betas, using_ages

def get_length_sides(betas, ages):
    min_age = min(ages)
    max_age = max(ages)
    vicinity = 5
    length = []
    print("INFO: Getting length of the sides...")
    for age in min_age, max_age:
        neighboring_betas = data.neighboring_betas(age, ages, betas, vicinity)
        sigma = stats.stdev(neighboring_betas)
        avg = stats.math_expectation(neighboring_betas)
        high_betas = avg + 2 * sigma
        low_betas = avg - 2 * sigma
        length.append(high_betas - low_betas)
    print("INFO: Getting length of the sides...DONE")
    return length