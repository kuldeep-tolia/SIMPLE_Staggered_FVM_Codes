! Finite-volume Navier Stokes solver
! SIMPLE algorithm on a staggered grid
! -------------------------------------------------------------------------------------------------------------------------------------------------
! Laminar Channel flow Problem:
! Velocity boundary conditions: no-slip stationary walls, uniform velocity (governed by Re) imposed at inlet boundary,
!                               zero normal velocity-gradient at outlet boundary
! Pressure boundary condition: Neumann boundary condition at all boundaries (requires a reference pressure value)
! Must ensure mass conservation at inlet and outlet boundaries
! -------------------------------------------------------------------------------------------------------------------------------------------------
! The pressure correction loop is solved using GS-SOR iterative method
! The diffusion fluxes are discretized using linear-profile assumption
! The convective fluxes are discretized using 1st order upwind scheme
! -------------------------------------------------------------------------------------------------------------------------------------------------
! To compile and run the program, use the following command at the terminal:
! > gfortran file_name.f90 -o executable_name.out
! > ./executable_name.out
! A .csv output file is generated to view the results in Paraview
! -------------------------------------------------------------------------------------------------------------------------------------------------

MODULE SOLVER_MEMORY
  IMPLICIT NONE

  INTEGER, PARAMETER :: rk = KIND(1.0D+00)
  INTEGER            :: nx, ny  
  INTEGER            :: u_iter_max, v_iter_max, pressure_iter_max
  INTEGER            :: clock_rate, clock_max, start_t, end_t

  DOUBLE PRECISION   :: global_tol, global_res
  DOUBLE PRECISION   :: dx, dy, xmax, xmin, ymax, ymin, xcoord, ycoord
  DOUBLE PRECISION   :: rho, mu, Re, in_velocity, in_flux
  DOUBLE PRECISION   :: alpha_u, alpha_v, alpha_p, mass_error

  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: u, v, p, u_old, v_old, pprime
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: u_intp, v_intp, p_intp
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: ap_u, ap_v, ap_p, ae, aw, an, as, mass_imb

  REAL(KIND = rk)    :: wall_time

END MODULE SOLVER_MEMORY

MODULE LDC_ROUTINES
  USE SOLVER_MEMORY

CONTAINS

  SUBROUTINE ALLOCATE_MEMORY()
    IMPLICIT NONE

    ALLOCATE (u(0:nx+1, 0:ny+1)     , v(0:nx+1, 0:ny+1)     , p(0:nx+1, 0:ny+1)     )
    ALLOCATE (u_old(0:nx+1, 0:ny+1) , v_old(0:nx+1, 0:ny+1) , pprime(0:nx+1, 0:ny+1))
    ALLOCATE (u_intp(0:nx+1, 0:ny+1), v_intp(0:nx+1, 0:ny+1), p_intp(0:nx+1, 0:ny+1))
    ALLOCATE (ap_u(0:nx+1, 0:ny+1)  , ap_v(0:nx+1, 0:ny+1)  , ap_p(0:nx+1, 0:ny+1)  )
    ALLOCATE (ae(0:nx+1, 0:ny+1)    , aw(0:nx+1, 0:ny+1)    , an(0:nx+1, 0:ny+1)    ) 
    ALLOCATE (as(0:nx+1, 0:ny+1)    , mass_imb(0:nx+1, 0:ny+1)                      )

  END SUBROUTINE ALLOCATE_MEMORY

  SUBROUTINE DEALLOCATE_MEMORY()
    IMPLICIT NONE

    DEALLOCATE (u     , v     , p     )
    DEALLOCATE (u_old , v_old , pprime)
    DEALLOCATE (u_intp, v_intp, p_intp)
    DEALLOCATE (ap_u  , ap_v  , ap_p  )
    DEALLOCATE (ae    , aw    , an    )
    DEALLOCATE (as    , mass_imb      )

  END SUBROUTINE DEALLOCATE_MEMORY

  SUBROUTINE SET_PARAMETERS()
    IMPLICIT NONE

    ! set number of grid points
    nx = 500                    ! set number of points in x-direction
    ny = 50                     ! set number of points in y-direction

    ! set domain size
    xmin = 0.0d0
    xmax = 10.0d0
    ymin = 0.0d0
    ymax = 1.0d0

    ! calculate grid spacing
    dx = (xmax - xmin) / FLOAT(nx)
    dy = (ymax - ymin) / FLOAT(ny)

    ! set governing parameters
    Re           = 100.0d0                                   ! set Reynolds number
    rho          = 1.0d0                                     ! set density of fluid
    in_velocity  = 1.0d0                                     ! set inflow velocity
    mu           = (rho * in_velocity * (ymax - ymin)) / Re  ! calculating dynamic viscosity
    in_flux      = in_velocity * rho * (ymax - ymin)         ! calculating inflow mass rate
    
    ! set relaxation factors
    alpha_u = 0.650d0             ! u-momentum relaxation factor
    alpha_v = 0.650d0             ! v-momentum relaxation factor
    alpha_p = 0.60d0              ! pressure correction relaxation factor

    ! set number of iterations
    global_tol = 1.0d-6            ! tolerance to converge SIMPLE algorithm
    u_iter_max = 30                ! max iterations of u-momentum solver
    v_iter_max = 30                ! max iterations of v-momentum solver
    pressure_iter_max = 250        ! max iterations of pressure correction solver
    
  END SUBROUTINE SET_PARAMETERS
  
  SUBROUTINE INITIALIZE_FVMSOLVER()
    IMPLICIT NONE

    ! initialize primary variables
    u     = 0.0d0
    v     = 0.0d0
    p     = 0.0d0

    ! set bc - north boundary
    u(:, ny+1) = 0.0d0
    v(:, ny+1) = 0.0d0

    ! set bc - south boundary
    u(:, 0)    = 0.0d0
    v(:, 1)    = 0.0d0

    ! set bc - west boundary
    u(1, :)    = in_velocity
    v(0, :)    = 0.0d0

    ! set bc - east boundary
    u(nx+1, :) = u(1, :)
    v(nx+1, :) = 0.0d0

    ! initialize all the solver coefficients
    ap_u = 1.0d0
    ap_v = 1.0d0
    ap_p = 1.0d0
    ae   = 0.0d0
    aw   = 0.0d0
    an   = 0.0d0
    as   = 0.0d0

    u_old = u
    v_old = v
    
  END SUBROUTINE INITIALIZE_FVMSOLVER

  SUBROUTINE U_MOMENTUM_SOLVER()
    IMPLICIT NONE

    INTEGER          :: i, j, k
    DOUBLE PRECISION :: de, dw, dn, ds, out_flux
    DOUBLE PRECISION :: me, mw, mn, ms, u_residual

    global_res = 0.0d0

    ! set coefficients for u-momentum discrete equation
    DO j = 1, ny
       DO i = 2, nx

          de         = mu * dy / dx
          dw         = mu * dy / dx
          dn         = mu * dx / dy
          ds         = mu * dx / dy
          
          me         = rho * dy * 0.50d0 * (u_old(i+1, j) + u_old(i, j)    )
          mw         = rho * dy * 0.50d0 * (u_old(i-1, j) + u_old(i, j)    )
          mn         = rho * dx * 0.50d0 * (v_old(i, j+1) + v_old(i-1, j+1))
          ms         = rho * dx * 0.50d0 * (v_old(i, j)   + v_old(i-1, j)  )

          IF (j .EQ. ny) dn = 2.0d0 * dn ! update coefficients along north boundary
          IF (j .EQ. 1 ) ds = 2.0d0 * ds ! update coefficients along south boundary

          ae(i, j)   = de + MAX(-me, 0.0d0)
          aw(i, j)   = dw + MAX( mw, 0.0d0)
          an(i, j)   = dn + MAX(-mn, 0.0d0)
          as(i, j)   = ds + MAX( ms, 0.0d0)

          ap_u(i, j) = ae(i, j) + aw(i, j) + an(i, j) + as(i, j) + &
               &       me - mw + mn - ms

       END DO
    END DO
    
    ap_u = ap_u / alpha_u
    
    ! iterating u-momentum equation
    DO k = 1, u_iter_max
       DO j = 1, ny
          DO i = 2, nx

             u(i, j) = ((1.0d0 - alpha_u) * u_old(i, j)) + ((1.0d0 / ap_u(i, j)) * &
                  &    (ae(i, j) * u(i+1, j) + &
                  &     aw(i, j) * u(i-1, j) + &
                  &     an(i, j) * u(i, j+1) + &
                  &     as(i, j) * u(i, j-1) + &
                  &     dy * (p(i-1, j) - p(i, j))))

          END DO
       END DO

       u(nx+1, :) = u(nx, :)    ! zero gradient at outlet
       
    END DO

    ! satisfy mass conservation
    out_flux = 0.0d0
    DO j = 1, ny

       out_flux = out_flux + rho * dy * u(nx+1, j)

    END DO

    u(nx+1, :) = in_flux / out_flux * u(nx+1, :) ! correction applied at outflow mass rate

    ! compute u-momentum residual
    u_residual = 0.0d0
    DO j = 1, ny
       DO i = 2, nx

          u_residual = u_residual + (u(i, j) - u_old(i, j))**2.0

       END DO
    END DO

    u_residual = SQRT(u_residual)

    global_res = u_residual

    WRITE(*, '(A, F15.8)') "u-momentum error = ", u_residual

  END SUBROUTINE U_MOMENTUM_SOLVER

  SUBROUTINE V_MOMENTUM_SOLVER()
    IMPLICIT NONE

    INTEGER          :: i, j, k
    DOUBLE PRECISION :: de, dw, dn, ds
    DOUBLE PRECISION :: me, mw, mn, ms, v_residual

    ! set coefficients for v-momentum discrete equation
    DO j = 2, ny
       DO i = 1, nx

          de         = mu * dy / dx
          dw         = mu * dy / dx
          dn         = mu * dx / dy
          ds         = mu * dx / dy

          me         = rho * dy * 0.50d0 * (u_old(i+1, j) + u_old(i+1, j-1))
          mw         = rho * dy * 0.50d0 * (u_old(i, j-1) + u_old(i, j)    )
          mn         = rho * dx * 0.50d0 * (v_old(i, j+1) + v_old(i, j)    )
          ms         = rho * dx * 0.50d0 * (v_old(i, j-1) + v_old(i, j)    )

          IF (i .EQ. nx) de = 2.0d0 * de ! update coefficients along east boundary
          IF (i .EQ. 1 ) dw = 2.0d0 * dw ! update coefficients along west boundary

          ae(i, j)   = (mu * dy / dx) + MAX(-me, 0.0d0)
          aw(i, j)   = (mu * dy / dx) + MAX( mw, 0.0d0)
          an(i, j)   = (mu * dx / dy) + MAX(-mn, 0.0d0)
          as(i, j)   = (mu * dx / dy) + MAX( ms, 0.0d0)

          ap_v(i, j) = ae(i, j) + aw(i, j) + an(i, j) + as(i, j) + &
               &       me - mw + mn - ms

       END DO
    END DO
    
    ap_v = ap_v / alpha_v

    ! iterating v-momentum equation
    DO k = 1, v_iter_max
       DO j = 2, ny
          DO i = 1, nx

             v(i, j) = ((1.0d0 - alpha_v) * v_old(i, j)) + ((1.0d0 / ap_v(i, j)) * &
                  &    (ae(i, j) * v(i+1, j) + &
                  &     aw(i, j) * v(i-1, j) + &
                  &     an(i, j) * v(i, j+1) + &
                  &     as(i, j) * v(i, j-1) + &
                  &     dx * (p(i, j-1) - p(i, j))))

          END DO
       END DO
    END DO

    ! compute v-momentum residual
    v_residual = 0.0d0
    DO j = 2, ny
       DO i = 1, nx

          v_residual = v_residual + (v(i, j) - v_old(i, j))**2.0

       END DO
    END DO

    v_residual = SQRT(v_residual)

    global_res = global_res + v_residual

    WRITE(*, '(A, F15.8)') "v-momentum error = ", v_residual

  END SUBROUTINE V_MOMENTUM_SOLVER

  SUBROUTINE PRESSURE_CORRECTION_SOLVER()
    IMPLICIT NONE

    INTEGER :: i, j, k

    ! set coefficients for pressure correction discrete equation
    DO j = 1, ny
       DO i = 1, nx

          ae(i, j) = rho * dy**2.0 / ap_u(i+1, j  )
          aw(i, j) = rho * dy**2.0 / ap_u(i  , j  )
          an(i, j) = rho * dx**2.0 / ap_v(i  , j+1)
          as(i, j) = rho * dx**2.0 / ap_v(i  , j  )

       END DO
    END DO

    ! update coefficients along boundaries
    ae(nx, :)  = 0.0d0
    aw(1,  :)  = 0.0d0
    an(:, ny)  = 0.0d0
    as(:, 1 )  = 0.0d0

    ap_p       = ae + aw + an + as
    ap_p(1, 1) = 1.0d40         ! setting reference cell pressure to initialized pressure

    ! initializing the pressure corrections to zero
    pprime     = 0.0d0

    ! computing mass imbalance in each cell before solving for pressure correction equation
    mass_imb = 0.0d0
    DO j = 1, ny
       DO i = 1, nx

          mass_imb(i, j) = rho * dy * (u(i+1, j) - u(i, j)) + &
               &           rho * dx * (v(i, j+1) - v(i, j))

       END DO
    END DO

    mass_error = SUM(mass_imb**2.0)
    mass_error = SQRT(mass_error)
    WRITE(*, '(A, F15.8)') "Mass imbalance before solving for pressure correction = ", mass_error

    ! solving for pressure correction equation
    DO k = 1, pressure_iter_max
       DO j = 1, ny
          DO i = 1, nx
             
             pprime(i, j) = pprime(i, j) + (1.50d0 / ap_p(i, j)) * &
                  &         (ae(i, j) * pprime(i+1, j) + &
                  &          aw(i, j) * pprime(i-1, j) + &
                  &          an(i, j) * pprime(i, j+1) + &
                  &          as(i, j) * pprime(i, j-1) - &
                  &          mass_imb(i, j)            - &
                  &          ap_p(i, j) * pprime(i, j))             

          END DO
       END DO
    END DO
    
    ! apply pressure corrections to pressure; must-be under-relaxed
    DO j = 1, ny
       DO i = 1, nx

          p(i, j) = p(i, j) + alpha_p * pprime(i, j)

       END DO
    END DO

    ! apply velocity corrections to u-velocity
    DO j = 1, ny
       DO i = 2, nx
          
          u(i, j) = u(i, j) + ((dy / ap_u(i, j)) * (pprime(i-1, j) - pprime(i, j)))

       END DO
    END DO
    
    ! apply velocity corrections to v-velocity
    DO j = 2, ny
       DO i = 1, nx
          
          v(i, j) = v(i, j) + ((dx / ap_v(i, j)) * (pprime(i, j-1) - pprime(i, j)))

       END DO
    END DO
    
    ! computing mass imbalance in each cell after solving for pressure correction equation
    DO j = 1, ny
       DO i = 1, nx

          mass_imb(i, j) = rho * dy * (u(i+1, j) - u(i, j)) + &
               &           rho * dx * (v(i, j+1) - v(i, j))

       END DO
    END DO
    mass_error = SUM(mass_imb**2)
    mass_error = SQRT(mass_error)
    
    WRITE(*, '(A, F15.8)') "Mass imbalance after solving for pressure correction = ", mass_error

    global_res = global_res + mass_error

    ! update velocity field
    u_old = u
    v_old = v
    
  END SUBROUTINE PRESSURE_CORRECTION_SOLVER

  SUBROUTINE POST_PROCESS_DATA()
    IMPLICIT NONE

    INTEGER :: i, j

    OPEN(UNIT = 3, FILE = "Converged_Results_LDC.csv")

    ! staggered grid has solution variables placed at different locations
    ! interpolating all interior solution variables at one location
    DO j = 2, ny
       DO i = 2, nx
          
          u_intp(i, j) = 0.50d0 * (u(i, j-1) + u(i, j))
          v_intp(i, j) = 0.50d0 * (v(i-1, j) + v(i, j))
          p_intp(i, j) = 0.250d0 * (p(i-1, j-1) + p(i-1, j) + p(i, j) + p(i, j-1))

       END DO
    END DO

    ! interpolating at boundaries
    u_intp(1, 2:ny)    = 0.50d0 * (u(1, 1:ny-1)    + u(1, 2:ny)   )
    u_intp(nx+1, 2:ny) = 0.50d0 * (u(nx+1, 1:ny-1) + u(nx+1, 2:ny))
    u_intp(2:nx, 1)    = u(2:nx, 0)
    u_intp(2:nx, ny+1) = u(2:nx, ny+1)
    
    v_intp(1, 2:ny)    = v(0, 2:ny)
    v_intp(nx+1, 2:ny) = v(nx+1, 2:ny)
    v_intp(2:nx, 1)    = 0.50d0 * (v(1:nx-1, 1)    + v(2:nx, 1)   )
    v_intp(2:nx, ny+1) = 0.50d0 * (v(1:nx-1, ny+1) + v(2:nx, ny+1))

    p_intp(1, 2:ny)    = 0.50d0 * (p(1, 1:ny-1)  + p(1, 2:ny) )
    p_intp(nx+1, 2:ny) = 0.50d0 * (p(nx, 1:ny-1) + p(nx, 2:ny))
    p_intp(2:nx, 1)    = 0.50d0 * (p(1:nx-1, 1)  + p(2:nx, 1) )
    p_intp(2:nx, ny+1) = 0.50d0 * (p(1:nx-1, ny) + p(2:nx, ny))

    ! interpolating at corners
    u_intp(1, 1)       = 0.0d0
    u_intp(1, ny+1)    = 0.0d0
    u_intp(nx+1, 1)    = 0.0d0
    u_intp(nx+1, ny+1) = 0.0d0
    
    v_intp(1, 1)       = 0.0d0
    v_intp(1, ny+1)    = 0.0d0
    v_intp(nx+1, 1)    = 0.0d0
    v_intp(nx+1, ny+1) = 0.0d0
    
    p_intp(1, 1)       = p_intp(2, 2)
    p_intp(1, ny+1)    = p(1, ny)
    p_intp(nx+1, 1)    = p(nx, 1)
    p_intp(nx+1, ny+1) = p(nx, ny)

    WRITE(3, *) "xcoord, ycoord, zcoord, p, u, v, w"
    
    DO i = 1, nx+1

       xcoord = dx * FLOAT(i-1)

       DO j = 1, ny+1

          ycoord = dy * FLOAT(j-1)
          WRITE(3, 11) xcoord, ycoord, 0.0, p_intp(i, j), u_intp(i, j), v_intp(i, j), 0.0
11        FORMAT (F8.5, ",", F8.5, ",", F8.5, ",", F8.5, ",", F8.5, ",", F8.5, ",", F8.5)

       END DO
    END DO

    CLOSE(UNIT = 3)          
    
  END SUBROUTINE POST_PROCESS_DATA
    
END MODULE LDC_ROUTINES

PROGRAM MAIN
  USE LDC_ROUTINES
  IMPLICIT NONE

  INTEGER :: iter

  CALL SYSTEM_CLOCK(start_t, clock_rate, clock_max)

  CALL SET_PARAMETERS()
  CALL ALLOCATE_MEMORY()
  CALL INITIALIZE_FVMSOLVER()

  global_res = 1.0d0
  iter = 1

  DO WHILE(global_res .GT. global_tol)

     WRITE(*, '(A, I6)') "SIMPLE Iteration = ", iter
     CALL U_MOMENTUM_SOLVER()
     CALL V_MOMENTUM_SOLVER()
     CALL PRESSURE_CORRECTION_SOLVER()
     WRITE(*, '(A, F15.8)') "Global residual = ", global_res
     WRITE(*, '(A)') "----------------------------------------------------------------------"
     iter = iter + 1
     
  END DO

  CALL POST_PROCESS_DATA()
  CALL DEALLOCATE_MEMORY()

  CALL SYSTEM_CLOCK(end_t)
  wall_time = REAL(end_t - start_t, KIND = rk) / REAL(clock_rate, KIND = rk)

  WRITE(*, '(A, F15.8)') "Governing parameter: Re = ", Re
  WRITE(*, '(A, F15.8)') "Governing parameter: mu = ", mu
  WRITE(*, '(A, F15.8)') "Governing parameter: rho = ", rho
  WRITE(*, '(A, F15.8)') "Governing parameter: U_in = ", in_velocity
  WRITE(*, '(A, I5, A, I5)') "The grid resolution used = ", nx, " x ", ny
  WRITE(*, '(A, F15.8, A)') "Program running time = ", wall_time, " seconds"
  WRITE(*, '(A)') "The SIMPLE iterations are computed. EXITING!!!"
  
END PROGRAM MAIN
