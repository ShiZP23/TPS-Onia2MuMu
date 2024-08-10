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
3. Import trigger results and save for possible trigger matching.
4. Harvest muons from tracks. Loop over muon pairs.
5. For muon pair candidataes, apply a crude selection with "opposite-charge criterion" and mass window cut. ( For $J/\psi$, consider mass range $[1.0, 4.0]$. For $\Upsilon$, consider mass range $[8.0, 12.0]$. )
6. Apply kinematic fitting for each pair to vertices. 
7. Store valid muon pairs by category. ($J/\psi$ or $\Upsilon$)
8. Use non-overlapping muon pairs to form $J/\psi$ and $\Upsilon$ candidates. (Careful not to directly store iterators! Use pointers instead.)
9. Store the candidates and corresponding muon pairs.

### Efficiency

### Systematics