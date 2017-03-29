# CS 205 Final Project: Mapping out Trajectories of Charged Defects

## Introduction
The promise of the nitrogen-vacany (NV) center in diamond as a system for implementing memory storage for a quantum computer has spawned interest in defect centers in related materials such as SiC. SiC is particularly interesting as it is a polymorphic material, exhibiting about 250 known polytypes, which imbues it with a degree of freedom unavailable in diamond. The three most common polytypes, 4H- and 6H-SiC and 3C-SiC, all have spin relaxation times ranging from 8 to 24 ms at 20 K (with 4H-SiC being the highest) and coherence persists up to room temperature [1]. In addition to the long spin coherence, a key feature is the ability to optically address (write in and read out) the spin states. However, much of the luminescence or emission of the defects is diverted into transitions involving scattering processes (and is not purely from the desired spin transitions) at ambient temperatures. Indeed, only about 4% of luminescence is from the desired transitions [2]. A potential solution is to place the defects near cavities on resonance with the desired transitions [3]. Positioning defects, however, is a non-trivial endeavor and we propose that the negatively charged silicon vacancy defect in 4H-SiC could be positioned through the probabilistic process of diffusion.

## Background
Work has already been done on obtaining barriers to diffusion and the vibrational frequencies that provide the directionality and time scale for the motion [4,5], but these studies lack a comprehensive understanding of the potential diffusion pathways and as a result do not perform the kinetic Monte Carlo simulations necessary to truly establish diffusion pathways and probabilities.

## Method
There are two main stages to the calculation of the trajectory probability maps. In the first stage, density functional theory calculations will be carried out to obtain the barriers to diffusion for the various pathways and in the second stage kinetic Monte Carlo simulations will be carried out to map out trajectories. We have existing code for the second part in Matlab, which will be converted to C code. The kinetic Monte Carlo simulations will be carried out both without and with Coulomb interaction and for 1D and 2D random walk test cases (to be compared with theoretical calculations). Without the coulomb interaction the calculation of the trajectories of multiple defects will be embarrassingly parallel.


1. A. L. Falk, B. B. Buckley, G. Calusine, W. F. Koehl, V. V. Dobrovitski, A. Politi, C. A. Zorman, P. X. L. Feng, and D. D. Awschalom, “Polytype control of spin qubits in silicon carbide,” Nature Communications, vol. 4, p. 1819, 05 2013.

2. I. Aharonovich, S. Castelletto, D. A. Simpson, C.-H. Su, A. D. Greentree, and S. Prawer, “Diamond-based single-photon emitters,” Reports on Progress in Physics, vol. 74, no. 7, p. 076501, 2011.

3. D. O. Bracher, X. Zhang, and E. L. Hu, “Selective purcell enhancement of two closely linked zero-phonon transitions of a silicon carbide color center,” arXiv:1609.03918, 2016.

4. E. Rauls, T. Frauenheim, A. Gali, and P. Deák, “Theoretical study of vacancy diffusion and vacancy-assisted clustering of antisites in sic,” Phys. Rev. B, vol. 68, p. 155208, Oct 2003.

5. X. Wang, M. Zhao, H. Bu, H. Zhang, X. He, and A. Wang, “Formation and annealing behaviors of qubit centers in 4h-sic from first principles,” Journal of Applied Physics, vol. 114, no. 19, p. 194305, 2013.
