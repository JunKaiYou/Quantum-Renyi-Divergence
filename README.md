# Quantum-Renyi-Divergences

This repository provides the codes used for numerical experiments on minimizing quantum Renyi divergences (partial results can be found in [this work](https://arxiv.org/abs/2109.06054)). In a nutshell, we propose two Polyak-type step sizes for mirror descent and apply the proposed methods to compute four quantum information quantities: Petz-Augustin information, sandwiched Augustin information, conditional sandwiched Renyi entropy, and sandwiched Renyi information.

### Contents
The contents of this repository are as follows:
- `Petz_Augustin_Information`: this file contains methods to compute Petz-Augustin information.
- `Sandwiched_Augustin_Information`: this file contains methods to compute sandwiched Augustin information.
- `Conditional_Sandwiched_Renyi_Entropy`: this file contains methods to compute conditional sandwiched Renyi entropy.
- `Sandwiched_Renyi_Information`: this file contains methods to compute sandwiched Renyi information.


### Requirements
- Install the MATLAB toolbox [QETLAB](http://www.qetlab.com/Main_Page), and unzip the file in your MATLAB scripts directory
- Install [Matrix Function Toolbox (MFT)](https://www.maths.manchester.ac.uk/~higham/mftoolbox/) 
- Install the MATLAB package given in [here](https://www.mathworks.com/matlabcentral/fileexchange/41621-fractional-matrix-powers-with-frechet-derivatives-and-condition-number-estimate) to compute the Frechet derivative of the power function[<sup>1</sup>](#refer-anchor-1)

<div id="refer-anchor-1"></div>
[1] Nicholas J. Higham and Lijing Lin. An improved Schur–Pade algorithm for fractional powers of a matrix and their frechet derivatives. SIAM J. Matrix Anal. Appl., 34(3):1341–1360, 2013.
