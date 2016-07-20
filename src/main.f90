program dmft_ed_zt
    use mpi
    use fdf
    use alloc, only: alloc_report, alloc_default
    use timer, only: t1_run, t2_run, timestamp, print_elapsed_time
    use dmft_params, only: dmft_params_read, &
                           nspin, norb, nbath, nsite, &
                           nloop, nev, nstep, &
                           nsector, sectors, nw, tbham, &
                           em_present, ek_present, vmk_present, &
                           em_init, ek_init, vmk_init, &
                           U, Up, Jex, Jp, Mu, scf_tol
    use dmft_grid, only: dmft_grid_init
    use dmft_lattice, only: dmft_lattice_init

    implicit none

    integer :: level
    double precision :: threshold
    character(len=100) :: input_file

    call mpi_setup

    if (iargc()==0) then
        input_file = "input.fdf"
    else
        call getarg(1,input_file)
    endif

    if (master) then
        write(6,"(a)") repeat("=",80)
        write(6,*)
        write(6,*) "Multi-orbital DMFT"
        write(6,"(a,I2)") "Number of processors = ", nprocs
        write(6,"(1x,A,A)") "Input file = ", trim(input_file)
        write(6,*)
        call timestamp("Start of run")
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

    ! start run timing
    t1_run = mpi_wtime(mpierr)
    
    ! Read run parameters from input fdf file
    call dmft_params_read
    call dmft_grid_init
    call dmft_lattice_init


    ! end run timing
    t2_run = mpi_wtime(mpierr)

    call alloc_report( printNow=.true. )

    if (master) then
        write(6,*)
        write(6,"(a)") repeat("=",80)
        write(6,*)
        call print_elapsed_time("Total Running Time", t1_run, t2_run)
        write(6,*)
        call timestamp("End of run")
        write(6,*)
        write(6,"(a)") repeat("=",80)
    endif

    call fdf_shutdown
    call mpi_shutdown
end program dmft_ed_zt
