#=

    FIGURE 2
    ========

    Patients 26 and 99
        (following Strobl et al. 2021)

    Show fits and predicted propotion of resistant cells for two patients.

    Fits for all other patients in the posterior are in the supplement.

=#

using RobustAdaptiveTherapy
using Plots

include("../defaults.jl")

fig2a = plot_patient_fit(26;n=1000,show_prediction=true)
fig2b = plot_patient_fit(99;n=1000,show_prediction=true)
fig2c = plot_patient_fit(36;n=1000,show_prediction=true)
fig2d = plot_patient_fit(85;n=1000,show_prediction=true)

fig2 = plot(fig2a, fig2b, fig2c, fig2d, widen=true,
    size=(1000,450), layout=(1,4))

savefig(fig2, split("$(@__FILE__).svg",".")[1] * ".svg")
fig2