# Geometric Local Parameterization Method for Generalized Hele–Shaw Free Boundary Problems with Source Terms
This repository contains the implementation of the methods described in our paper (https://arxiv.org/abs/2607.02880). The project focuses on solving generalized Hele-Shaw problems with source terms on manifolds which are identified by point clouds. 

The provided codes display the boundary evolution dynamics and reproduce the main numerical examples from the paper:

(1) FBP_exact_source.m - Circular case with exact radial source term shown in Figure 2

(2) FBP_radial_approx.m - Circular case with approximate radial source term shown in Figure 4

(3) FBP_tumor.m - Applications to tumor growth model shown in Figure 6 and Figure 7 with different values for the parameters, D_1 and D2.

(4) FBP_nonradial_approx.m - Hele-Shaw free boundary problems with approximate nonradial source term shown in Figure 9.

If you use this code in your work, please cite:

@article{zhang2025geometric, title={A Geometric Local Parameterization Method for Generalized Hele-Shaw Free Boundary Problems with Source Terms}, author={Zhang, Zengyan and Hao, Wenrui and Harlim, John},
journal={arXiv preprint arXiv:2607.02880}, year={2026} }
