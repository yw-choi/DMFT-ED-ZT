module dmft
!==============================================================================
! controller module for dmft calculation.
! this module takes care of the flow of the program,
! and the actual works are done by other modules called by this module.
!==============================================================================
    use mpi
    use dmft_grid
    use dmft_params
    use dmft_lattice
    use dmft_green, only: G, G0, Sigma, G_prev, &
                          dmft_green_init, &
                          update_local_green_ftn, &
                          update_weiss_ftn, &
                          dump_green, &
                          find_density

    use impurity_solver, only: impurity_solver_init, solve, solver_post_processing

    use utils, only: die
    use timer, only: t1_loop, t2_loop, print_elapsed_time, t1_run, t2_run
    use alloc, only: alloc_report, re_alloc, de_alloc, alloc_default
    use fdf

    implicit none

    public :: &
        dmft_init,            &
        dmft_loop,            &
        dmft_post_processing, &
        dmft_finalize          

    integer ::  &
        iloop     ! DMFT loop count

    logical :: &
        converged  ! is DMFT loop converged?
contains

    subroutine dmft_init
        integer :: level
        double precision :: threshold
        character(len=100) :: input_file

        call mpi_setup
        if (iargc()==0) then
            input_file = "input.fdf"
        else
            call getarg(1,input_file)
        endif
        t1_run = mpi_wtime(mpierr)

        if (master) then
            write(6,"(a)") repeat("=",80)
            write(6,*)
            write(6,*) "Multi-orbital DMFT"
            write(6,*) "Number of processors = ", nprocs
            write(6,"(1x,A,A)") "Input file = ", trim(input_file)
            write(6,*)
            write(6,"(a)") repeat("=",80)
        endif

        call fdf_init(input_file, 'fdf.out')

        level = fdf_get('alloc_report_level', 2)
        threshold = fdf_get('alloc_report_threshold', 0.d0)
        call alloc_report( level=level, file='alloc_report', &
                           threshold=threshold, &
                           printNow=.false. )

        call alloc_default(copy=.false., shrink=.true.)

        ! Read general DMFT parameters independent of solver.
        call dmft_params_read
        call dmft_grid_init
        call dmft_lattice_init
        call dmft_green_init

        call impurity_solver_init
    end subroutine dmft_init

    subroutine dmft_loop
        integer :: ispin, ia

        if (master) then
            write(6,"(a)") repeat("=",80)
            write(*,*)
            write(*,*) "Start of DMFT SCF loop"
            write(*,*)
            write(6,"(a)") repeat("=",80)
        endif

        converged = .false.
        iloop = 0
        do while(.not.converged.and.iloop.le.nloop)
            iloop = iloop + 1
            if (master) then
                write(6,"(a)") repeat("=",80)
                write(*,*)
                write(*,*) "DMFT SCF loop ", iloop
                write(*,*)
                write(6,"(a)") repeat("=",80)
            endif

            t1_loop = mpi_wtime(mpierr)

            do ia=1,na
                call solve(iloop,ia,G0(:,:,:,ia), Sigma(:,:,:,ia))
            enddo

            call update_local_green_ftn
            call update_weiss_ftn
            call find_density
            call dump_green 

            t2_loop = mpi_wtime(mpierr)
            call loop_end
        enddo

        if (master) then
            write(*,"(a)") repeat("=",80)
            write(*,*)
            write(*,*) "End of DMFT SCF loop"
            write(*,*)
            write(*,"(a)") repeat("=",80)
        endif
    end subroutine dmft_loop

    subroutine loop_end
        integer :: iw, ia, iorb
        double precision :: diff, diffsum

        ! Test DMFT convergence.
        diff = sum(abs(G_prev(:,:,:,:)-G(:,:,:,:)))
        diff = diff / (na*norb*nspin*nwloc)

        call mpi_allreduce(diff,                 &
                           diffsum,              &
                           1,                    &
                           mpi_double_precision, &
                           mpi_sum,              &
                           comm,                 &
                           mpierr)

        if (diffsum < scf_tol) then
            converged = .true.
        endif

        if (master) then
            write(*,*)
            write(*,"(1x,a,I4,a,E12.5)") "loop ",iloop, &
                                         " done. scf diff = ",diffsum
            call print_elapsed_time(" Elapsed time in a loop     ",&
                                    t1_loop,t2_loop)
            call print_elapsed_time(" Elapsed time from the start",&
                                    t1_run,t2_loop)
            if (converged) then
                write(*,*) "DMFT loop has converged within ",iloop," iterations."
            endif
            write(*,*)
        endif

    end subroutine loop_end

    subroutine dmft_post_processing
        integer :: iw, iorb, ispin, ia
        character(len=100) fn
        double precision :: zq

        if (master) then
            write(*,*) 
            write(*,*) "Quasiparticle weight"
            write(*,*) 
            write(*,*) "    Z(ia,ispin,iorb)"
            write(*,*) "    ---------------------------"
            do ia=1,na
                do ispin=1,nspin
                    do iorb=1,norb
                        zq = aimag(sigma(1,iorb,ispin,ia))&
                            /omega(1)                        
                        zq = 1.d0/(1.d0-zq)
                        write(*,"(a,I2,a,I4,a,I3,a,F12.6)") &
                            "     Z(",ia,",",ispin,",",iorb,",) = ",zq
                    enddo
                enddo
            enddo
            write(*,*) 
        endif

        call mpi_barrier(comm,mpierr)

        call solver_post_processing

    end subroutine dmft_post_processing

    subroutine dmft_finalize
        t2_run = mpi_wtime(mpierr)

        call alloc_report( printNow=.true. )

        if (master) then
            write(6,*)
            write(6,"(a)") repeat("=",80)
            write(6,*)
            call print_elapsed_time(" Total Running Time", t1_run, t2_run)
            write(6,*)
            write(6,*) "End of run."
            write(6,*)
            write(6,"(a)") repeat("=",80)
        endif

        call fdf_shutdown
        call mpi_shutdown
    end subroutine dmft_finalize

end module dmft
