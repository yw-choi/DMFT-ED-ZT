module impurity_solver
!==============================================================================
! Wrapper module for impurity solvers. 
! This module separates the specific implementation of the impurity solver with
! DMFT module. 
!==============================================================================
    use fdf
    use utils
    use dmft_params
    use dmft_grid
    use ed_solver
    use mpi

    implicit none

    character(len=100) :: solver

contains

    subroutine impurity_solver_init
        if (master) then
            write(*,*) "Initializing impurity solver..."
        endif

        solver = fdf_get("DMFT.Solver", "ED")

        select case(solver)
            case ("ED")
                call ed_init
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select
    end subroutine impurity_solver_init

    subroutine solve(iloop,ia,G0,Sigma)
        integer, intent(in) :: iloop, ia
        double complex, intent(in) :: G0(nwloc,norb,nspin)
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)

        if (master) then
            write(*,"(1x,a,I1,a,I2,a)") "impurity problem for ia = ",ia
        endif

        select case(solver)
            case ("ED")
                call ed_solve(iloop,ia,G0,Sigma)
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select

    end subroutine solve

    subroutine solver_post_processing
        select case(solver)
            case ("ED")
                call ed_post_processing
            case default
                call die("impurity_solver", "Specfieid solver is not implemented.")
        end select

    end subroutine solver_post_processing

end module impurity_solver
