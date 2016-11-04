# NIST Industrial Wireless Analysis
This repository contains code for analyzing industrial wireless sounder measurements
and the generation of wireless scenarios.

![EM Wave](images/EMwave.jpg)

**DISCLAIMER**: Certain commercial equipment, instruments, or materials are identified in this paper in order to specify the experimental procedure adequately. Such identification is not intended to imply recommendation or endorsement by NIST, nor is it intended to imply that the materials or equipment identified are necessarily the best available for the purpose.

## Repository Structure
### Main components
* meas_analysis:    Includes files for analyzing RF Sounder measurements of various factories. 
  Data may be found at http://doi.org/10.18434/T44S3N
  * The main M code entry point is *estimate_channel_cwd.m* for processing of NIST channel impulse response MAT files
  * The file *stats2rfnestdp.m* is used to reduce average impulse responses to usable size within our testbed

* tap_reduction:    Library for generating an N-sampled impulse response for RF channel emulator integration

### Supporting capabilities
* www_scripts:      Python scripts used to generate file listings used by the www data landing page

## Meas_Analysis M-Code Library
### Main analysis code (M code)
MATLAB computer code used for the analysis of RF propagation measurement data.  Analysis includes the following:
* Gain versus Distance
* Line-of-sight versus Non-Line-of-Sight
* Rician K-factor
* First andnSecond Order delay spread
* Figure production for publication
* Statistical approximations
* Tabulation of results

### +reporting module (M code)
Code used for the generation of the NIST Technical Report XYZW latex-based report. Reporting code includes the following:
* Figure production for publication
* Statistical approximations
* Tabulation of results into Latex tables
* Automation of full latex-based reporting




