# SIM_Reconstruction

## Description

Structured illumination microscopy (SIM) reconstruction algorithm with computational sectioning.  More details can be found in the following paper:

<a href="https://ieeexplore.ieee.org/document/8579117" target="_blank">Computational Super-Sectioning for Single-Slice Structured Illumination Microscopy.</a>, <br />
IEEE Transactions on Computational Imaging, vol. 5 no. 2, pp. 240-250, 2019.  <br />
E. Soubies, and M. Unser.

## Requirements

The code requires the GlobalBioIm library v1.1.2 (or more recent releases) <br />
https://biomedical-imaging-group.github.io/GlobalBioIm/

## Repository content

The repository is organized as follows.

- The script **SimScript2D.m** contains the main code of the method proposed in [1].
- The script **SimuSIM2D.m** allow to generate 2D-SIM data (used in ScriptSimuExample.m).
- The script **ScriptFig3.m** reproduces the Figure 3 of [1].
- The script **ScriptSimuExample.m** is an example on a small simulated example without out-of-focus signal.
- The script **ScriptPatternFromFairSIM.m** provides an example on how to generate images of the patterns given the xml file that can be obtained with <a href="https://www.fairsim.org/" target="_blank"> FairSIM </a>
- Folder **Data** contains raw-SIM data with corresponding PSF and patterns generated as explained in ScriptPatternFromFairSIM.m
- Folder **Utils** contains auxilliary functions

