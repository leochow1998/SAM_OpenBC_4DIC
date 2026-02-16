# SAM-OpenBC-4DIC: System for Atmospheric Modeling, with open boundary condition and 4D initialization
Author: Leo Chow - leonardochow1998@link.cuhk.edu.hk <br/>
Date: 2026-02-16 <br/>
Version: 6.10.4.1

This is a modified version of System for Atmospheric Modeling version 6.10.4 (SAM) to suit the experimental need in Chow et al. (2026). The modifications are:
1. `dampling_lateral.f90`: Meridional lateral sponge layer nudging zonal wind and temperature close to the meridional boundaries to a zonal averaged reference value, meridional wind, vertical wind, and microphysics variables to zero. The method follows Klemp and Lilly (1978).
2. `open_boundary_y.f90`: Meridional open boundary condition. The implementation follows Klemp and Wilhelmson (1978). The code modifications are documented in `/SRC/Change_openybc.txt`.
3. `read4dcaminput.f90`: Read pre-processed input nc files from CAM to initialize the SAM model, instead of initialize it from a spatially uniform vertical profile or previous SAM runs. 

**Caution**: This version of SAM is only tested with the configuration stated in the paper. Therefore, the author do not gurantee the code will run in other configurations and settings. 

### Data:
Namelist files are stored in `Qobs10_SAM.000`. <br/>
Since the data used for initializing SAM in the paper is too large in size, please send an email to Leo Chow (leonardochow1998@link.cuhk.edu.hk) if you want to obtain a copy of it.

### To cite this repository:
Paper: Chow et al. (2026) (Chow, T. N., C. Y. Tam, and T. S. E. Chung, 2026: Investigating Convective Influences on Tropical Cyclone Genesis in a Super-parameterized General Circulation Model Aquaplanet Experiment. *Journal of Atmosphertic Science*, (in preperation). <br/>
Software: [![DOI](https://zenodo.org/badge/1155858269.svg)](https://doi.org/10.5281/zenodo.18653589) Chow et al. (2026) (Chow, T. N., C. Y. Tam, and T. S. E. Chung, 2026: SAM-OpenBC-4DIC: System for Atmospheric Modeling, with open boundary condition and 4D initialization (Version 6.10.4.1). Zenodo, [http://doi.org/10.5281/zenodo.18653590](http://doi.org/10.5281/zenodo.18653590)
