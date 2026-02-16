# SAM-OpenBC-4DIC: System for Atmospheric Modeling with Open Boundary Conditions and 4D Initialization
Author: Leo Chow - leonardochow1998@link.cuhk.edu.hk <br/>
Date: 2026-02-16 <br/>
Version: 6.10.4.1

This is a modified version of the System for Atmospheric Modeling version 6.10.4 (SAM) from Khairoutdinov and Randall. (2003), customized for the experimental requirements in Chow et al. (2026). Key modifications include:
1. `damping_lateral.f90`: Implements a meridional lateral sponge layer. This nudges zonal wind and temperature toward a zonally averaged reference value near the meridional boundaries, while meridional wind, vertical wind, and microphysics variables are nudged toward zero. The methodology follows Klemp and Lilly (1978).
2. `open_boundary_y.f90`: Implements meridional open boundary conditions following Klemp and Wilhelmson (1978). Detailed code modifications are documented in `/SRC/Change_openybc.txt`.
3. `read4dcaminput.f90`: Allows the model to read pre-processed NetCDF input files from CAM for initialization, rather than initializing from a spatially uniform vertical profile or previous SAM runs.

Comments within the source code provide further implementation details.

**Caution**: This version of SAM has only been tested with the configurations described in the associated paper. The author does not guarantee that the code will function correctly in other configurations or settings.

### Data:
Namelist files are stored in `Qobs10_SAM.000`. <br/>
Due to the large size of the initialization data used in the paper, please contact Leo Chow (leonardochow1998@link.cuhk.edu.hk) to obtain a copy.

### To cite this repository:
Paper: Chow et al. (2026) (Chow, T. N., C. Y. Tam, and T. S. E. Chung, 2026: Investigating Convective Influences on Tropical Cyclone Genesis in a Super-parameterized General Circulation Model Aquaplanet Experiment. *Journal of Atmosphertic Science*, (in preperation).) <br/>
Software: [![DOI](https://zenodo.org/badge/1155858269.svg)](https://doi.org/10.5281/zenodo.18653589) Chow et al. (2026) (Chow, T. N., C. Y. Tam, and T. S. E. Chung, 2026: SAM-OpenBC-4DIC:  System for Atmospheric Modeling with Open Boundary Conditions and 4D Initialization (Version 6.10.4.1). Zenodo, [http://doi.org/10.5281/zenodo.18653590](http://doi.org/10.5281/zenodo.18653590).)
