Problem Description:  

-> Flow inside a channel/duct is another canonical problem.  
-> Here, a FVM-based Navier-Stokes solver is presented utilizing SIMPLE algorithm on a staggered grid.  
-> Boundary conditions:  
- Velocity boundary conditions: no-slip stationary walls, uniform velocity (governed by Re) imposed at inlet boundary, zero normal velocity-gradient at outlet boundary   
- Pressure boundary condition: Homogenous Neumann boundary condition at all boundaries $\frac{\partial p}{\partial n} = 0$  

-> The pressure-correction loop is solved using Gauss-Seidel SOR iterative method.  
-> The diffusive fluxes are discretized using a linear-profile assumption.  
-> The convective fluxes are discretized using $1^{st}-$ order upwind scheme.  
-> The results for the $Re = 100$ case are presented here.  
