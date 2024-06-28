# High-Resolution-Spectroscopy-Analysis-of-HD-149026-b
This repository contains coding scripts used in Rafi+ (2024) (accepted to AJ), which are used to analyze high-resolution spectroscopy (HRS) data of the hot Saturn HD 149026 b with CARMENES/NIR.

The HRS data is taken from the publicly available database of CARMENES that can be accessed from [here](http://caha.sdc.cab.inta-csic.es/calto/jsp/searchform.jsp). For convenience, the near-infrared (NIR) data (in FITS files) from fiber A, which is used in the analysis, is provided in this repository under `data/spectrum/` folder. This is for **transmission spectroscopy** data only.

## Dependencies

The code was successfully run using Python v3.9.18, whilst the older nor newer version has not been tested yet.

The code requires several Python packages to be installed. These are:
- `astropy`
- `scipy`
- `PyAstronomy`
- `pandas`
- `tqdm`
- `h5py`
- `batman`
- `petitRADTRANS` (pRT2 version; pRT3 has not been tested yet)
- `pyfastchem`

Please refer to the documentation of each package for the installation procedure, but they should be able to be installed by simply using `pip` command:
```
pip install <package>
```

Specifically for `petitRADTRANS` pRT2, after installation, the opacity data needs to be downloaded (provided in the [documentation](https://petitradtrans.readthedocs.io/en/latest/content/installation.html)) and then add an environment variable to the downloaded opacity data folder:
```
import os
os.environ['pRT_input_data_path'] = 'relative/path/to/input_data'
```

Please see the [documentation](https://petitradtrans.readthedocs.io/en/latest/content/installation.html) for details. 

For this analysis, I used the H2O opacity from HITEMP 2010 and POKAZATEL line list ([Polyansky+ 2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.2597P/abstract)) computed using HELIOS-K ([Grimm+ 2021](https://ui.adsabs.harvard.edu/abs/2021ApJS..253...30G/abstract)), as well as for HCN opacity already provided in the petitRADTRANS documentation. The computed transmission spectrum templates (specifically for HD 149026 b's assumed atmospheric condition) are provided under the `models/` folder. It is up to the users if they want to compute models on their own using their own opacity tables. Regardless, here I provide the necessary materials to recreate the analysis in Rafi+ (in prep.).

## Steps of Main Analysis

The analysis is divided into four main steps:

### Data reduction

The corresponding script for this step is `data_reduction_original.ipynb`. 

The code consists of data selection (masking low S/N frames and spectral orders), sky emission removal with sky emission templates from [Oliva+ (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...581A..47O/abstract) and [Rousselot+ (2000)](https://ui.adsabs.harvard.edu/abs/2000A%26A...354.1134R/abstract), residual blaze correction, and outliers removal. In addition, a wavelength stability check is also performed. However, this case, in particular, does not need any additional wavelength correction since the wavelength solution is already stable, hence the code does not consist of an algorithm to correct the wavelength solution.

The code also uses the `batman` package [(Kreidberg 2015)](https://ui.adsabs.harvard.edu/abs/2015PASP..127.1161K/abstract) to build the light curve model for the corresponding transit of this planet, which will act as a weight for each frame when we cross-correlate them with the model.

Finally, to remove the dominating stellar and telluric lines, infamously known for being major obstacles in HRS studies for exoplanet atmospheres, the code uses the SysRem algorithm described in [Tamuz+ (2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.356.1466T/abstract) which will fit and correct any systematics present in the data. After passing the data into SysRem, the 'clean' data is obtained and ready to be cross-correlated.

To compute the systemic velocity $V_{sys}$ of the system, the radial velocity (RV) value of the star is also calculated. The circular radial velocity equation is given by:
$$RV_s(\phi)=K_s\sin(2\pi\phi)+V_{sys}$$
Knowing $RV_s$ (stellar RV), $K_s$ (stellar orbital velocity), and $\phi$ (orbital phase), $V_{sys}$ can then be fitted to the radial velocity curve ($RV$ as a function of time or orbital phase $\phi$).

### Model generation

The corresponding script for this step is `model_generation.ipynb`. 

The code uses `pyfastchem` to calculate chemical equilibrium abundances and `petitRADTRANS` to calculate the model spectrum (single and multiple trace species). In addition, the code can also calculate simplified constant species abundances typically used when prior knowledge of the atmospheric abundances is not known. Although `petitRADTRANS` can compute both transmission and emission spectrum, this code is written only for transmission.

In calculating the transmission spectrum, the code will first use `petitRADTRANS` to compute the spectrum, then broaden it to both the FWHM of the instrument resolution and the planetary rotational motion (assuming circular orbit with $v\sin{i}=2\pi r/P$). Both broadening functions are given in `function.py`. When normalization is required, the code will use high-pass (also given in `function.py`) and Gaussian filter to compute the pseudo-continuum, then either divide or subtract the continuum from the model (depending on the user's request).

### Cross-correlation
The corresponding script for this step is `cross_correlation.ipynb`.

Typical high-resolution spectroscopy studies for exoplanet atmospheres use the cross-correlation technique, where the SysRem residuals are cross-correlated with the model spectrum, to boost and search the planet's atmospheric signal. The code will first cross-correlate the data and model along an $RV_p$ grid, then build the so-called S/N map by interpolating the resulting CCF values into a grid of $K_p$ (planetary orbital velocity) and $V_{rest}$ values, where the latter (a variable to account for any additional systematic velocity after accounting the orbital and systemic velocities) follows:
$$RV_p(\phi)=K_p\sin(2\pi\phi)+V_{sys}+V_{rest}+V_{bary}$$
where subscript $p$ is for the planet. $V_{bary}$ is the barycentric velocity during the observation. 

### Pre-processing

The corresponding script for this step is `pre_processing.ipynb`.

To perform retrieval analysis, where an exact match between data and model is favorable, pre-processing the model is therefore necessary. This is to account for the fact that SysRem typically alters the planet signal, resulting in a planetary spectrum that deviates significantly from the model ([Gibson+ 2022](https://academic.oup.com/mnras/article/512/3/4618/6510825)). Pre-processing the model aims to mimic this SysRem effect in the model spectrum. After pre-processing, the script will perform retrieval analysis using a simple grid search to find the best-fit orbital and systemic velocity (which in this case is $V_{rest}$). To reproduce the figure comparing likelihood between including and excluding frame #44 (Figure 14 in Rafi+ 2024), use `plot_pre_processing.ipynb` script.

## Additional Analysis

Several other analyses to confirm whether the signal is real (in the case where signal S/N is low, such as this one) are also performed. These are the Welch-t test and injection test. The former can be found in `welch-t.ipynb` whilst the latter can be found in `cross_correlation.ipynb`. Please see Rafi+ (2024) for details about these tests. Additionally, to check whether telluric lines may somehow contribute to the observed (water) signal, cross-correlation between the SysRem residuals and telluric model template is also performed. The script is given in `cross_correlation_telluric.ipynb`.

On the other hand, another test is performed to see how different reduction parameter values, particularly $\eta$ (it sets how wide the line wings around strong telluric lines will be masked), can affect the resulting signal S/N. This test can be found in `data_reduction_original_eta.ipynb` and `cross_correlation_eta.ipynb`.

In Rafi+ (2024), we tried to separate the transit into two halves to check whether the red-shifted $V_{rest}$ of the H2O signal that we found originated solely from the first half of the transit where H2O flows from the night-side to the day-side along the planet's leading terminator. The script for this analysis can be found in `sanity_check_frame.ipynb`.

## How to Use

The code is a just-shift-and-enter-code, meaning the user can simply run the code without changing or adjusting anything, with the following sequence: 
1. Run `data_reduction_original.ipynb`. It will export the clean data (after SysRem) to `data/hdf5/`.
2. Run `model_generation.ipynb`. It will export the model spectra to `models/`. This step can be skipped if the users want to use the already computed models provided in the `models/` folder.
3. Run `cross_correlation.ipynb` to see the resulting S/N maps.
4. Run `pre_processing.ipynb` to retrieve and constrain the orbital and systemic velocity values.

Indeed, the parameter values (e.g., in the data reduction, model generation, or cross-correlation) are tunable and they can be changed by the users.

Please feel free to contact me (Rafi) via email at salirafi8@gmail.com if you have any questions or post an issue in the repo!
