Problem Description:  

-> Lid-driven cavity is a standard test case with reference data available in the literature.  
-> Here, a FVM-based Navier-Stokes solver is presented utilizing SIMPLE algorithm on a staggered grid.  
-> Boundary conditions:  
- Velocity bc: no-slip stationary walls $\left( u = v = 0 \right)$, top plate imposed with uniform velocity $\left( u = 1, v = 0 \right)$    
- Pressure boundary condition: Homogenous Neumann boundary condition at all walls $\frac{\partial p}{\partial n} = 0$  

-> The pressure-correction loop is solved using Jacobi iterative method, which is parallelized using OpenMP.  
-> The diffusive fluxes are discretized using a linear-profile assumption.  
-> The convective fluxes are discretized using a hybrid scheme of $1^{st}-$ order upwind scheme and CDS.  
-> The results are compared with the data available in the literature.  

Reference: Ghia, U. K. N. G., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. Journal of computational physics, 48(3), 387-411.  
