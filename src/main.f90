program dmft_ed_zt
    use mpi
    use fdf
    use timer, only: t1_run, t2_run, timestamp, print_elapsed_time

    implicit none

    character(len=100) :: input_file

    call mpi_setup

    if (iargc()==0) then
        ! default input fdf file
        input_file = "input.fdf"
    else
        ! or read file name from the first command line argument
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

    ! start run timing
    t1_run = mpi_wtime(mpierr)

    ! ==============================================================
    ! RUN DMFT SUBROUTINE
    ! ==============================================================
    call dmft
    ! ==============================================================

    ! end run timing
    t2_run = mpi_wtime(mpierr)

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
