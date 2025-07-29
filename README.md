# Anisotropy-analysis-JBC-W2025

### Written by C. Milano 
Author's contact Info: cmilano@uiowa.edu Additional PoC: frederick-stern@uiowa.edu
Last Updated: July 29th, 2025

## User Agreement and Approval for Usage
Users of this software are not allowed re-distribute the files of this repository (including all subdirectories) anywhere. Furthermore, usage of this software is only granted with the expressed approval of IIHR - Hydroscience & Engineering. Please contact Dr. Frederick Stern (frederick-stern@uiowa.edu) for approval.

## Description:
A Python code for extracting turbulence characteristic quantities at two vortex cores. This program is designed to output anisotropic Reynolds stresses, eigenvalues of the Reynolds stress tensor and anisotropic Reynolds stress tensor, turbulence invariants, and Reynolds stress ellipsoid.

## Notes:
This program may require the use of basic Python packages, such as NumPy, Pandas, SciPy, etc. These should be easily installed using the Pip package manager. The same requirements.txt used in 'Vortex Core Turbulence Analysis' should cover all the necessary packages.

<br>

This program does not require additional inputs from the command line, but some files are required to input the Reynolds stress tensor and the mean velocity vector.

<br>

The output of the programs are written into the same directory as the code location. 


## Examples

The simplest example is to have 4 files formatted like the files shown in *Example_data*. The files in *Example_data* can be replaced with the actual values of EFD/CFD mantaining the same format to obtain the desired outpts. The program can be used with the following command:
```
python anisotropy_analysis_abv1_abv2.py
```
<br>

## Theoretical background
The code uses the Reynolds stress tensor:

<img width="210" height="74" alt="image" src="https://github.com/user-attachments/assets/c61ee411-8cbb-4c4d-bfe5-b6c844e65a59" />

and the mean velocity vector:

<img width="162" height="37" alt="image" src="https://github.com/user-attachments/assets/747f64c9-38de-4b3b-b704-308e857b137b" />

to calculate turbulent related quantities. The anisotropic Reynolds stress tensor is calculated as:

<img width="191" height="106" alt="image" src="https://github.com/user-attachments/assets/be522dfe-90ea-4b19-a32e-3f75884e7b14" />

The eigenvalues of the Reynolds stress tensor are calculates using: 

<img width="193" height="48" alt="image" src="https://github.com/user-attachments/assets/f9139669-81a7-4204-ab94-e170d9adff6b" />

and ordered from largest to smallest. The corresponding eigenvalues of the Anisotropic Reynolds stress tensor are:

<img width="170" height="42" alt="image" src="https://github.com/user-attachments/assets/ad13cefa-533e-4fb8-86db-abcc6dbbee66" />

The turbulence invariants are calculated as:

<img width="270" height="61" alt="image" src="https://github.com/user-attachments/assets/958e9525-6ccc-47c5-a75e-407c16eb16aa" />

The Reynolds stress ellipsoid is defined by the expected value equation:

<img width="101" height="29" alt="image" src="https://github.com/user-attachments/assets/de963d44-f6dc-40cb-af33-ba0a9043ba5b" />

or given constant 𝛼 ≥ 0 with solution given as:

<img width="282" height="50" alt="image" src="https://github.com/user-attachments/assets/e31d334c-e307-4a8e-94f6-1ceb84880951" />

This equation represents an ellipsoid centered in the mean velocity with the half-lengths of the principal axed being 𝛼 times the square root of the Reynolds stress tensor eigenvalues. Assuming 𝛼=1 gives the surface of the ellipsoid. 




