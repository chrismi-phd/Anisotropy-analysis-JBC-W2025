# Anisotropy-analysis-JBC-W2025

### Written by C. Milano 
 Contact Info: cmilano@uiowa.edu
 
### July 29, 2025

## Description:
A Python code for extracting turbulence characteristic quantities at two vortex cores. This program is designed to output anisotropic Reynolds stresses, eigenvalues of the Reynolds stress tensor and anisotropic Reynolds stress tensor, turbulence invariants, and Reynolds stress ellipsoid.

## Notes:
This program may require the use of basic Python packages, such as NumPy, Pandas, SciPy, etc. These should be easily installed using the Pip package manager. The same requirements.txt used in 'Vortex Core Turbulence Analysis' should cover all the necessary packages.

<br>

This program does not require additional inputs from the command line, but some files are required to input the Reynolds stress tensor and the mean velocity vector.

<br>

The output of the programs are written into the same directory as the code location. 


## Examples

The simplest example is to have 4 files formatted like the files shown in *Example_data*. The files in *Example_data* can be replaced with the actual values of EFD/CFD mantaining the same format to obtain the desired outpts. The program can be  ran with the following command:
```
python anisotropy_analysis_abv1_abv2.py
```

tau_ij = mu * (∂u_i/∂x_j + ∂u_j/∂x_i)
<br>

