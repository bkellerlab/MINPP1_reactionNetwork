# MINPP1_reactionNetwork

this repository contains supporting information to the manuscript:<br>
**"Stable isotopomers of myo-inositol to uncover the complex MINPP1-dependent inositol phosphate network"**, <br>
M. Nguyen Trung, S. Kieninger, Z. Fandi, D. Qiu, G. Liu, A. Saiardi, H. Jessen, B. G. Keller, D. Fiedler <br>
preprint: https://doi.org/10.1101/2022.08.29.505671 (2022)
<br>

In our work, we numerically analyzed experimentally determined progress curves for the MINPP1 mediated dephosphorylation of InsP<sub>5</sub>[2OH] and InsP<sub>6</sub>.
This repository contains the corresponding Python3 scripts that were used to extract the kinetic rates from the experimental data sets. 


<br>
<br>
<h2>Overview </h2>

<h4>1D model system </h4>

All analysis was performed with Python3. The files "fitted_exp_data.py" contain all fit functions that where used to fit the scaled experimental data.
The files "minimization_routine.py"contain the code for the entire minimization routine including initial guess and reconstruction of the progress curves
from the evaluated rates.
The directory tree is as follows:

```
Code
  |
  ├── InsP5[2OH]
  |       |   fitted_exp_data.py
  |       |   minimization_routine.py
  |
  |      
  └── InsP6
          |   fitted_exp_data.py
          |   minimization_routine.py
```
