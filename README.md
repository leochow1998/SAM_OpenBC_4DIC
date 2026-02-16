# SAM_OpenBC_4DIC - System for Atmospheric Modeling, with open boundary condition and 4D initialization
Author: Leo Chow - leonardochow1998@link.cuhk.edu.hk

This is a modified version of System for Atmospheric Modeling version 6.10.4 (SAM) to suit the experimental need in Chow et al. (2026). The modifications are:
1. dampling_lateral.f90: Meridional lateral sponge layer nudging zonal wind and temperature close to the meridional boundaries to a zonal averaged reference value, meridional wind, vertical wind, and microphysics variables to zero. The method follows Klemp and Lilly (1978)
2. open_boundary_y.f90: Meridional open boundary condition. The implementation follows Klemp and Wilhelmson (1978). The code modifications are documented in /SRC/Change_openybc.txt
3. read4dcaminput.f90: Read pre-processed input nc files from CAM to initialize the SAM model, instead of initialize it from a spatially uniform vertical profile or previous SAM runs. 

Caution: This version of SAM is only tested with the configuration stated in the paper. Therefore, the author do not gurantee the code will run in other configurations and settings. 

Data used to initialize the SAM experiment in Chow et al. (2026): 

To cite this repository:
Paper: Chow et al. (2026) (Chow, T. N., C. Y. Tam, and T. S. E. Chung, 2026: Investigating Convective Influences on Tropical Cyclone Genesis in a Super-parameterized General Circulation Model Aquaplanet Experiment. Journal of Atmosphertic Science (in preperation)
Software: 
