# 2025-PGforLQRwithStochastic
This repository contains the code from our paper "Convergence Guarantees of Model-free Policy Gradient Methods for LQR with Stochastic Data"

For Figure 1: please run Figure1.m first and then use FigureGeneration.m to plot the figure.

For Figure 2: run Figure2.m

For Figure 3: run RunPGD_3_withoutVR.m, RunPGD_5_withoutVR.m, RunPGD_3_withVR.m, RunPGD_5_withVR.m, RunPGD_3_withoutVR_ex.m. Run ProcessData.m to process the data generated for the black lines (remove useless part due to divergence). Finally run FigureGeneration.m

For Figure 4 and Figure 5: run NPG3_ada.m NPG5_ada.m NPG3_nonada.m NPG5_nonada.m. Then run ProcessData.m to process data to remove NaN entries. Finally run Figure4Generation.m or Figure5Generation.m to generate the figures.
