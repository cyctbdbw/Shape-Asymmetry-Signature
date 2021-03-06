Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

# Shape-Asymmetry-Signature
![image](https://user-images.githubusercontent.com/48049448/169631022-87e6cedf-8588-4226-9f9f-eff81ee1476c.png)

This repository provides scripts for the analysis in "The individuality of shape asymmetries of the human cerebral cortex". 
See https://www.biorxiv.org/content/10.1101/2021.09.14.460242v2 for details.

Python code will be provided soon.

For the heritability analysis, also see
Arnatkeviciute. A. et al., Genetic influences on hub connectivity of the human connectome. Nat Commun 12, 4237 (2021).
https://www.nature.com/articles/s41467-021-24306-2

Dependances: 
1. gifti toolbox (https://github.com/gllmflndn/gifti.git): this toolbox is used when generating functional connectivity
2. PALM (Permutation Analysis of Linear Models; https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM): this toolbox keep subjects' family structures intac when shuffling the subjects.
3. Permutation Inference for Canonical Correlation Analysis (https://github.com/andersonwinkler/PermCCA)
4. toolbox of correcting p-values for multiple comparison (https://au.mathworks.com/matlabcentral/fileexchange/61659-multicmp) 
5. plot the surface ROI boundary(https://github.com/StuartJO/plotSurfaceROIBoundary)
6. Shape-DNA (https://github.com/Deep-MI/BrainPrint-legacy/blob/master/fs_shapeDNA.py) and see the relevant papers below:
(1) M. Reuter, F.-E. Wolter, N. Peinecke, Laplace–Beltrami spectra as ‘Shape-DNA’ of surfaces and solids. Computer-Aided Design 38, 342-366 (2006). https://doi.org/10.1016/j.cad.2005.10.011
(2) M. Reuter, F. E. Wolter, M. Shenton, M. Niethammer, Laplace-Beltrami Eigenvalues and Topological Features of Eigenfunctions for Statistical Shape Analysis. Comput Aided Des 41, 739-755 (2009). https://doi.org/10.1016/j.cad.2009.02.007
