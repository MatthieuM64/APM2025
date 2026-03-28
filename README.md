# Metastability of liquid state of the active Potts model

Codes used in the scientific publication: S. Chatterjee, M. Karmakar, M. Mangeat, H. Rieger, and R. Paul, <i>Stability of discrete-symmetry flocks: Sandwich state, traveling domains and motility-induced pinning</i>, <a href='https://journals.aps.org/pre/abstract/10.1103/1r19-ryx9'>Phys. Rev. E <b>112</b>, 064115 (2025)</a>. A preprint is available on <a href='https://arxiv.org/abs/2507.08187'>arXiv</a>.</br></br>

A C++ code to compute the numerical simulations of the microscopic model (with/without the insertion of the droplet) and a C++ code to compute the numerical solutions of the hydrodynamic equations (with the insertion of the droplet) are available in this repository. A Python code for each c++ code is also available to generate movies of the system dynamics. </br></br>

<b>Exportations:</b> Density and state snapshots (microscopic simulations), density and magnetization snapshots (hydrodynamic solutions).</br>
<b>Compile:</b> g++ filename.cpp -fopenmp -lgsl -lgslcblas -lm -O3 -s -o filename.out.</br>
<b>Run:</b> ./filename.out -parameter=value.</br>
<b>Generate the movie:</b> python filename.py -parameter=value.</br>

## Droplet insertion (microscopic simulations)

<b>List of parameters (<i>APM_micro_droplet_omp.cpp</i>):</b> beta, D, epsilon, rho0, Rd, rhod, LX, LY, tmax, init, ran, threads (details as comments in the code).</br>
<b>List of parameters (<i>figure_APM_micro_droplet.py</i>):</b> beta, D, epsilon, rho0, Rd, rhod, LX, LY, tmax, init, ran, NCPU, movie (details as comments in the code).</br>

## Metastability of the liquid under small diffusion (microscopic simulations)

<b>List of parameters (<i>APM_micro_omp.cpp</i>):</b> beta, D, epsilon, rho0, LX, LY,tmax, init, ran, threads (details as comments in the code).</br>
<b>List of parameters (<i>figure_APM_micro.py</i>):</b> beta, D, epsilon, rho0, LX, LY, tmax, init, ran, NCPU, movie (details as comments in the code).</br>

## Droplet insertion (hydrodynamic solutions)

<b>List of parameters (<i>APM_hydro_omp.cpp</i>):</b> beta, epsilon, rho0, Rd, rhod, LX, LY, dt, tmax, dx, init, threads (details as comments in the code).</br>
<b>List of parameters (<i>figure_APM_hydro.py</i>):</b> beta, epsilon, rho0, Rd, rhod, LX, LY, tmax, init, NCPU, movie (details as comments in the code).</br>
