subroutine dmft

    use mpi
    use timer, only: t1_loop, t2_loop, timestamp
    use dmft_grid, only: dmft_grid_init, nwloc, omega
    use eigpair, only: eigpair_t
    use impurity_hamiltonian, only: initial_h_imp_params
    use dmft_lattice, only: dmft_lattice_init
    use dmft_green, only: cluster_hybridization_ftn, &
                          cluster_green_ftn, &
                          local_green_ftn
    use ed_projection, only: fit_h_imp_params 
    use dmft_params, only: dmft_params_read, &
                           nspin, norb, nbath, nloop

    implicit none

    ! data variables
    double complex, allocatable :: &
        G_cl(:,:,:),     & ! Gimp(nwloc,norb,nspin) 
        G_loc(:,:,:),     & ! Gloc(nwloc,norb,nspin)
        D_cl(:,:,:),       & ! D_cl(nwloc,norb,nspin)
        G0(:,:,:),        & ! G0(nwloc,norb,nspin) Weiss field
        Sigma(:,:,:)        ! Sigma(nwloc,norb,nspin) self energy

    double precision, allocatable :: &
        occ(:,:),        & ! occ(norb,nspin) occupancy from the local Green's ftn
        em(:,:),         & ! em(norb,2) impurity levels
        ek(:,:),         & ! ek(nbath,2) bath levels
        vmk(:,:,:)         ! vmk(norb,nbath,2) impurity-bath hybridization 

    ! ground state
    type(eigpair_t) :: gs

    ! local variables
    logical :: converged
    integer :: iloop
    double precision :: diff
    character(len=200) :: msg
    
    ! @TODO remove the following debugging variables
    integer :: iw
    
    ! Read run parameters from input fdf file
    call dmft_params_read

    ! matsubara frequency grid
    call dmft_grid_init

    ! tight-binding hamiltonian
    call dmft_lattice_init

    allocate(G_cl(nwloc,norb,nspin), G_loc(nwloc,norb,nspin))
    allocate(D_cl(nwloc,norb,nspin))
    allocate(G0(nwloc,norb,nspin), Sigma(nwloc,norb,nspin))
    allocate(occ(norb,nspin))
    allocate(em(norb,2), ek(nbath,2), vmk(norb,nbath,2))

    ! initial impurity hamiltonian parameters
    call initial_h_imp_params( em, ek, vmk )

    call cluster_hybridization_ftn( ek, vmk, D_cl )

    ! Main DMFT loop
    converged = .false.
    dmftloop: do iloop=1,nloop

        if (master) then
            write(6,"(a)") repeat("=",80)
            write(*,*)
            write(msg,"(A,I4)") "DMFT SCF loop ", iloop
            call timestamp(msg)
            write(*,*)
            write(6,"(a)") repeat("=",80)
        endif

        t1_loop = mpi_wtime(mpierr)

        call ground_state( em, ek, vmk, gs )

        call cluster_green_ftn( em, ek, vmk, gs, G_cl )

        call local_green_ftn( em, D_cl, G_cl, G_loc )
        open(unit=123,file="green.dump",status="replace")
        do iw=1,nwloc
            write(123,"(5F12.6)") omega(iw), real(G_cl(iw,1,1)), aimag(G_cl(iw,1,1)),&
                        real(G_loc(iw,1,1)), aimag(G_loc(iw,1,1))
        enddo
        close(123)

        call test_convergence( G_cl, G_loc, diff, converged )

        if (master) then
            write(*,*) "|G_cl - G_loc| = ", diff
            if (converged) then
                write(*,*) "DMFT loop has converged within ",iloop," iterations."
                exit dmftloop
            endif
        endif

        ! find new em, ek, vmk 
        call fit_h_imp_params( em, ek, vmk, D_cl, G_cl, G_loc )

        call cluster_hybridization_ftn( ek, vmk, D_cl )

        t2_loop = mpi_wtime(mpierr)
    enddo dmftloop

end subroutine dmft

subroutine test_convergence( G_cl, G_loc, diff, converged )
    use mpi
    use dmft_params, only: nspin, norb, nbath, nw, scf_tol
    use dmft_grid, only: nwloc

    double complex, intent(in) :: &
        G_cl(nwloc, norb, nspin), &
        G_loc(nwloc, norb, nspin)

    logical, intent(out) :: converged
                                 
    double precision :: diff_loc, diff

    diff_loc = sum(abs(G_cl-G_loc)) / (norb*nbath)
    call mpi_allreduce( diff_loc, diff, 1, mpi_double_precision, &
                        mpi_sum, comm, mpierr )

    converged = diff < scf_tol

end subroutine test_convergence
