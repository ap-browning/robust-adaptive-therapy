## Code associated with the preprint "Adaptive therapy under parametric, structural, and measurement uncertainty" available on bioRxiv

### Getting Started

Download this repository and activate the Julia environment using 
```
pkg> activate
```

The code required to produces the results is a Julia package RobustAdaptiveTherapy.
```
add /path/to/folder/RobustAdaptiveTherapy
```

You should then be able to call
```
using RobustAdaptiveTherapy
```
to access functions provided by the package. The first time you call this it will automatically download the data from the repository https://github.com/reneebrady/IADT_PCaSC into the `RobustAdaptiveTherapy/data` folder.

### Reproducing results

* **Fits:** The fits for individual patients are stored as `JLD2` files in the respective subfolder for each model within the `RobustAdaptiveTherapy/fits` folder. Therefore, all calls to `patient_fit(idx)` will by default load the saved fit unless the keyword argument `load_saved` is set to `false`.

* **Figures:** Code to produce the figures is available in the `figures` folder. For example, to reproduce Figure 2, run the file `figures/fig2/fig2.jl`.

### Supplementary fits

The respective `RobustAdaptiveTherapy/fits` subfolders also contain `.png` files showing the data and patient fit (essentially Figure 2) for patients with sufficient data available.