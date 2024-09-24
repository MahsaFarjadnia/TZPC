# Robust Data-Driven Tube-Based Zonotopic Predictive Control with Closed-Loop Guarantees
This repository contains all the necessary scripts and files to reproduce the results of Example 1 and Example 2 as presented in the following conference paper:

M. Farjadnia, A. Fontan, A. Alanwar, M. Molinari, and K. H. Johansson, "Robust Data-Driven Tube-Based Zonotopic Predictive Control
with Closed-Loop Guarantees," Accepted at 63rd Conference on Decision and Control, 2024, ([Link](https://arxiv.org/abs/2409.14366) to the arXiv version).

This work proposes a robust data-driven tube-based zonotopic predictive control (TZPC) approach for discrete-time linear systems, designed to ensure stability and recursive feasibility in the presence of bounded noise. The proposed approach consists of two phases. In an initial learning phase, we provide an over-approximation of all models consistent with past input and noisy state data using zonotope properties. Subsequently, in a control phase, we formulate an optimization problem, which, by integrating terminal ingredients is proven to be recursively feasible. Moreover, we prove that implementing this data-driven predictive control approach guarantees robust exponential stability of the closed-loop system.

## File Descriptions:

### Example 1:

#### "Example_1_TZPC.m":
This script employs the TZPC algorithm to compute optimal control inputs for Example 1 as provided in the referenced paper.

####  "Example_1_Plots.m":
This script generates the plots presented in the paper.
  
### Example 2:

#### "Example_2_State_Space_Model.m":
This script computes the discrete-time system matrices used in the paper. It also includes detailed explanations of the building model.

#### "Example_2_TZPC.m":
This script employs the TZPC algorithm to compute optimal control inputs for Example 2 as provided in the referenced paper.

####  "Example_2_Plots.m":
This script generates the plots presented in the paper.

### Prerequisites to run "Example_1_TZPC.m" and "Example_2_TZPC.m":
- Ensure the [MPT](https://www.mpt3.org/) and [Yalmip](https://yalmip.github.io/) toolbox are installed, along with the [MOSEK](https://www.mosek.com/products/academic-licenses/) solver.
- Add the 'Saved Workspace' folder to the MATLAB path.
- Add the 'RequiredFiles_TZPC' folder and its subfolders to the MATLAB path.

### Prerequisites to run "Example_1_Plots.m" and "Example_2_Plots.m":
- Add the 'Saved Workspace' folder to the MATLAB path.
- Add the 'RequiredFiles_TZPC' folder and its subfolders to the MATLAB path.
  
## Acknowledgments
Please note that portions of this code are derived from the [CORA](https://tumcps.github.io/CORA/) toolbox.

## Citation
The referenced paper's BibTeX is as follows:
```sh
@inproceedings{Farjadnia2024robust,
 title={A Robust Data-Driven Tube-Based Zonotopic Predictive Control with Closed-Loop Guarantees},
 author={Mahsa Farjadnia, Angela Fontan, Amr Alanwar, Marco Molinari, and Karl Henrik Johansson},
 booktitle={2024 IEEE 63rd Conference on Decision and Control (CDC)},
 pages={-},
 year={2024},
 address={Milan, Italy}
 }
```

