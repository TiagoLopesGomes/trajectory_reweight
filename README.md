# EnsembleLab BME Multiprocessing

Ensemble analysis and reweighting of protein trajectories using Bayesian/Maximum Entropy (BME) methods, with a focus on SAXS (Small Angle X-ray Scattering) data. It utilizes multiprocessing and trajectory files to speed up computations.

## Description

The `ensemblelab_bme_multiproc.py` script is designed to:

1. Load and analyze trajectories (currently GROMACS xtc files)
2. Calculate various structural metrics (Rg, end-to-end distance, Dmax, etc.)
3. Compute SAXS profiles for each frame of the trajectory
4. Perform BME reweighting using experimental SAXS data
5. Optimize the θ parameter for reweighting
6. Generate plots comparing original and reweighted ensembles

## Dependencies

- NumPy
- Pandas
- SciPy
- Scikit-learn
- MDAnalysis
- Matplotlib
- Seaborn
- tqdm
- kneed
- concurrent
- BME (Bayesian/Maximum Entropy reweighting)

Additional software:
- Pepsi-SAXS (for SAXS calculations)

Note:
- The BME.py and BME_tools.py are provided in the repository for importing BME functions 
and retrieved from https://github.com/KULL-Centre/BME

## Usage

1. Ensure all dependencies are installed.
2. Place your trajectory and topology files in a `./traj/` directory:
   - `./traj/{NAME}_trajectory.xtc`
   - `./traj/top_{NAME}.pdb`
3. Prepare your experimental SAXS data file.
4. Modify the `NAME` variable in the script to match your system.
5. Run the script:
```python ensemblelab_bme_multiproc.py```

## Key Functions

- `cm_dist`: Calculate distance between centers of mass of two selections
- `dmax`: Calculate maximum pairwise distance in a selection
- `eed`: Calculate end-to-end distance
- `calc_traj`: Calculate various metrics for each frame of the trajectory
- `kde`: Kernel Density Estimation for distribution plotting
- `autoblock`: Perform block analysis for error estimation
- `traj2saxs`: Calculate SAXS profile for each frame of the trajectory
- `iBME`: Perform iterative BME reweighting for a given θ value
- `find_optimal_theta`: Find the optimal θ value using the knee method

## Output

The script generates several plots:
1. χ² vs φ_eff plot for θ optimization
2. SAXS profiles (original vs. reweighted)
3. Kratky plots
4. Residuals plot
5. Distribution plots for Rg, end-to-end distance, Dmax, and center of mass distance

## Notes

- The script uses multiprocessing to speed up SAXS calculations and BME optimizations.
- Experimental SAXS data is preprocessed and errors are adjusted using BIFT (not included in this script).
- The optimal θ value is chosen using the knee method on the χ² vs φ_eff plot, or can be manually set.

## References

This script is adapted from this Jupyter notebook (incorporation of multiprocessing and using a trajectory to speed up analysis):
https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb

