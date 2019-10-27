import stats
import data

def get_low_high_betas(betas, ages):
    min_age = min(ages)
    max_age = max(ages)
    high_betas = []
    low_betas = []
    vicinity = 5
    if max_age - min_age <= vicinity:
        print("INFO: The range of ages is no larger than the given vicinity")
        sigma = stats.stdev(betas)
        avg = stats.math_expectation(betas)
        high_betas.append(avg + 2 * sigma)
        low_betas.append(avg - 2 * sigma)
        return low_betas, high_betas
    else:
        print("INFO: Getting upper and lower beta values...")
        for age in range(min_age, max_age+1, vicinity):
            neighboring_betas = data.neighboring_betas(age, ages, betas, vicinity)
            sigma = stats.stdev(neighboring_betas)
            avg = stats.math_expectation(neighboring_betas)
            high_betas.append(avg + 2 * sigma)
            low_betas.append(avg - 2 * sigma)
        print("INFO: Getting upper and lower beta values...DONE")
        return low_betas, high_betas