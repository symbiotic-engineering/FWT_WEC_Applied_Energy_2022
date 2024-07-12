# FWT-WEC-Cost-Dynamics
This repository is an open source codebase that optimizes the energy storage capacity, losses, and cost for a system of combined wave energy converters with offshore wind-wave farms with different renewable energy penetration levels. The floating wind turbine (FWT) used is **Siemens SWT-3.6 MW** , and a spar platform is used for mooring foundation. The spar buoy serves as a point-absorber wave energy converter (WEC).


**Context**

The project is part of research in the [Symbiotic Engineering Analysis (SEA) Lab](https://sea.mae.cornell.edu/)


**Citation**
- Kluger, J. M., Haji, M. N., & Slocum, A. H. (2023). The power balancing benefits of wave energy converters in offshore wind-wave farms with energy storage. Applied Energy, 331, 120389.


**Authors**
- Jocelyn M. Kluger, jociek@alum.mit.edu
- Maha Haji, maha@cornell.edu
- Alexander H. Slocum, slocum@mit.edu


**File Structure**
- `Dynamics`: Contains the main code `run_code` and other helper functions to simulate the dynamics of the WEC and FWT.
- `Optimization`: Contains the main code `check_results` and other helper functions to calculate and optimize the energy storage capacity, power supply and demand, losses, and cost.
- `Data`: Contains the data needed as input for main codes (already included in the main folders, Dynamics and Optimization). 
- `pubs`: Published manuscripts with the results from the codebase.


**How to use**
- To run the dynamics code, open and run `run_code.m`. Make sure all the helper functions and input data are in the same folder as the main code.
- To run the optimization code, open and run `check_results.m`. You can change the renewable energy penetration level as desired. Make sure all the helper functions and input data are in the same folder as the main code.


**Dependencies**
The following packages are used in this code:
- MATLAB
- Optimization Toolbox
- Statistics and Machine Learning Toolbox


**License**

This project is released open-source under the MIT License. The hydrodynamics data was generated from WAMIT and WEC-Sim, which is publicly available.
