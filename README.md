## General Info

Measuring the cross section of prompt $J/\psi+J/\psi+\Upsilon$ prodution in pp collision on CMS.

Aimed at providing measurement result for studying parton correlation in proton and background estimation for BSM searches.

## Development 

Using CMSSW_13_0_2. Inherited most of the code from AliceQuen/Onia2MuMu. (cr. Qin Junkai at Dept. of Physics, Tsinghua University)

### Contributors
* Wang Chi (Eric100911), undergraduate at Zhili College, Tsinghua University.
* Cheng Xing, undergraduate at Zhili College, Tsinghua University.
* Shi Zhenpeng, undergraduate at Dept. of Physics, Tsinghua University.

Supervised under Prof. Hu Zhen at Dept. of Physics, Tsinghua University.

## Overview

### Event Selection Procedure

1. Match $\mu^+\mu^-$ pairs using vertex fitting. (Track geometry only.)

2. Create $J/\psi$ and $\Upsilon$ candidates using $\mu^+\mu^-$ pairs. (Dynamics required. Mass window, pT selection and other restrictions required.)

3. Matching $J/\psi$ and $\Upsilon$ candidates from one single vertex . (Track geometry only. May check $c\tau$ distribution.)

### Efficiency

### Systematics

## Code Framework

### Event Selection Procedure

1. `void LoadMC()` Load MC results if `doMC == true`.
2. Some other code for initialization 

### Efficiency

### Systematics