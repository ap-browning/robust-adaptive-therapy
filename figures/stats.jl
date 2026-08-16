# Statistics quoted in the paper

using RobustAdaptiveTherapy
using DataFrames
using StatsBase

## Relapse 

    # Resamples to use
    n = 10000

    # Determine if relapse will eventuale
    relapse = q -> model_stable_equilibria(q)[2] > 0.0

    # Across all patients   
    Rall = count(relapse.(sample_joint_posterior(n))) / n

    # Across each of the patients in Figure 2
    idx = [26,99,36,85]
    Ridx = [count(relapse.(sample(patient_fit(i),n))) / n for i = idx]
    DataFrame(patient = idx, proportion = Ridx)

## Average duration

idx = patient_fit_list()

average_num_years = mean([maximum(patient_data(i)[1]) for i = idx]) / 365