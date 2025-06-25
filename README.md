# acf-galaxy-clustering

## Overview
This repository contains the latest version of the code used for the angular correlation function (ACF) analysis presented in:

* ["Galaxy clustering measurements out to redshift z∼8 from Hubble Legacy Fields"](https://ui.adsabs.harvard.edu/abs/2024MNRAS.528..898D/abstract)
* ["Galaxy clustering at cosmic dawn from JWST/NIRCam observations to redshift z∼11"](https://ui.adsabs.harvard.edu/abs/2024MNRAS.533.2391D/abstract)

The analysis measures galaxy clustering using the two-point ACF, incorporating both bootstrap and jackknife resampling techniques. This implementation includes an example application on [GOODS-S JADES v1](https://archive.stsci.edu/hlsp/jades) data, with random catalogs generated following the methodology described in the JWST projects above.

## Data Description

We omit the data reduction procedures as every case study requires different approaches. In this example, we provide:

- **Galaxy data**: Selected to be in the redshift range z = 5-6
- **Magnitude cut**: Absolute magnitude in UV < -17.0  
- **Survey field**: Observed in [GOODS-S JADES v1](https://archive.stsci.edu/hlsp/jades)
- **Quality control**: Cross-checked within the boundaries of the corresponding scientific, RMS, and segmentation maps
- **Random points**: Generated following the methodology described in ["Galaxy clustering at cosmic dawn from JWST/NIRCam observations to redshift z∼11"](https://ui.adsabs.harvard.edu/abs/2024MNRAS.533.2391D/abstract)

## Repository Structure

```
├── galaxy_acf_analysis.ipynb    # Example usage: main analysis notebook for measuring ACF 
├── ACF_calc.py                  # Core correlation function calculations
├── coordinate_data/             # Directory containing coordinate files
│   ├── coo_data_mag-17.0_z5.0.txt
│   └── coo_rand_mag-17.0_z5.0.txt
└── README.md                    # This file
```

## Requirements

### Python Dependencies
```
numpy
matplotlib
scipy
tqdm
```

### Custom Module
- `ACF_calc.py`: Contains bootstrap and jackknife correlation function estimators

## Methods

### Angular Correlation Function
The analysis uses the [Landy-Szalay estimator](https://ui.adsabs.harvard.edu/abs/1993ApJ...412...64L/abstract)

## Output

The analysis produces:

- Observed angular correlation function (ACF) measurements: ω(θ)
- Symmetric uncertainty estimates are derived from the covariance matrix, using error estimates from both resampling methods
- Visualization plots showing:
  - Survey field with data and random point distributions
  - ACF measurements with error bars
<p align="center">
  <img src="https://github.com/user-attachments/assets/7ccce6f2-247c-4649-b792-1511076b4b0e" alt="Survey Field" width="45%" />
  <img src="https://github.com/user-attachments/assets/e9fcd677-a9b3-474f-9172-8c37accaad01" alt="ACF Results" width="45%" />
</p>

## Citations

If you use this code, please cite:
- ["Galaxy clustering measurements out to redshift z∼8 from Hubble Legacy Fields"](https://ui.adsabs.harvard.edu/abs/2024MNRAS.528..898D/abstract)
- ["Galaxy clustering at cosmic dawn from JWST/NIRCam observations to redshift z∼11"](https://ui.adsabs.harvard.edu/abs/2024MNRAS.533.2391D/abstract)

## Contact
For questions about the code or methodology, please contact [ndalmasso](nicolo.dalmasso1@gmail.com).

**Note**: This is a demonstration example. For production analysis, ensure proper data validation, error handling, and parameter optimization for your specific dataset.
