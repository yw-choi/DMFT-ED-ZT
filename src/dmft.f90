subroutine dmft

    use mpi
    use timer, only: t1_loop, t2_loop, timestamp
    use dmft_grid, only: dmft_grid_init, nwloc
    use ed_basis, only: basis_t
    use impurity_hamiltonian, only: initial_h_imp_params
    use eigpair, only: eigpair_t
    use dmft_lattice, only: dmft_lattice_init
    use dmft_green, only: calc_Delta
    use dmft_params, only:dmft_params_read, &
                           nspin, norb, nbath, nsite, nloop, nsector, &
                           maxnstep, nw, tbham, em_present, ek_present, &
                           vmk_present, read_gf_from_file, &
                           U, Up, Jex, Jp, Mu, scf_tol, &
                           sectors

    implicit none

    ! data variables
    double complex, allocatable :: &
        G_old(:,:,:),     & ! G_prev(nwloc,norb,nspin) prev local Green's ftn
        G(:,:,:),         & ! G(nwloc,norb,nspin) current local Green's ftn
        Delta_old(:,:,:), & ! Delta_old(nwloc,norb,nspin) hybridization function
        Delta(:,:,:),     & ! Delta(nwloc,norb,nspin) hybridization function
        G0(:,:,:),        & ! G0(nwloc,norb,nspin) Weiss field
        Sigma(:,:,:)       ! Sigma(nwloc,norb,nspin) self energy

    double precision, allocatable :: &
        occ(:,:),        & ! occ(norb,nspin) occupancy from the local Green's ftn
        em(:,:),         & ! em(norb,2) impurity levels
        ek(:,:),         & ! ek(nbath,2) bath levels
        vmk(:,:,:)         ! vmk(norb,nbath,2) impurity-bath hybridization 

    ! ground state's basis, E0, gs 
    type(basis_t) :: basis
    
    type(eigpair_t) :: gs

    ! local variables
    logical :: &
        converged

    integer :: &
        iloop

    character(len=200) :: msg
    
    ! Read run parameters from input fdf file
    call dmft_params_read

    ! matsubara frequency grid
    call dmft_grid_init

    ! tight-binding hamiltonian
    call dmft_lattice_init

    allocate(G_old(nwloc,norb,nspin), G(nwloc,norb,nspin))
    allocate(G0(nwloc,norb,nspin), Sigma(nwloc,norb,nspin))
    allocate(Delta_old(nwloc,norb,nspin),Delta(nwloc,norb,nspin))
    allocate(occ(norb,nspin))
    allocate(em(norb,2), ek(nbath,2), vmk(norb,nbath,2))

    ! initial impurity hamiltonian parameters
    call initial_h_imp_params(em, ek, vmk)

    call calc_Delta(ek, vmk, Delta)

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

        call ground_state(em, ek, vmk, basis, gs)

        ! @TODO
        ! call cluster_green_ftn(em,ek,vmk, nev, eigpairs, G)

        ! @TODO test convergence with Green's function or Delta?
        ! call test_convergence(G_old, G, converged)

        if (converged) then
            if (master) then
                write(*,*) "DMFT loop has converged within ",iloop," iterations."
            endif
            exit dmftloop
        endif

        ! Delta_old = Delta
        ! calculate a new Delta from em, Delta_old, G
        ! call new_delta(em, Delta_old, G, Delta)

        ! @TODO
        ! find new em, ek, vmk by fitting Delta_cl to Delta
        ! call fit_delta_cl(Delta, em, ek, vmk)

        t2_loop = mpi_wtime(mpierr)
    enddo dmftloop

end subroutine dmft
