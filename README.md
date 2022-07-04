# Frenkel-Holstein-Approach-Steady-State-Spectra
Calculates aggregate steady state spectra based on Frenkel-Holstein Hamiltonian used extensively by Spano and co-workers (Chem Rev 2018 has all the necessary references). Codes are in MATLAB. The codes take into account finite temperature as well as nonzero disorder. For nonzero disorder, use at least 10000 configurations.

To reproduce the figures in the paper, copy over variables_chlbenz_abs.m to variables.m and run average_over_config.m in the MATLAB command line. Follow a similar step for emission in chlorobenzene and absorption in THF solvents.

The coherence part of the code is not final, ignore this section.

The original set of codes were writte by Divij Mishra and Pranav Kasetty (also authors on the above paper) and the codes are available at-
https://github.com/PranavKasetty/Linear-Agg-Temp-PL.git
