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
        G_cl(:,:,:),      & ! G_cl(nwloc,norb,nspin) 
        G_loc(:,:,:),     & ! G_loc(nwloc,norb,nspin)
        D_cl(:,:,:)         ! D_cl(nwloc,norb,nspin)

    double precision, allocatable :: &
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
    
    ! Read run parameters from input fdf file
    call dmft_params_read

    ! matsubara frequency grid
    call dmft_grid_init

    ! tight-binding hamiltonian
    call dmft_lattice_init

    allocate(G_cl(nwloc,norb,nspin), G_loc(nwloc,norb,nspin))
    allocate(D_cl(nwloc,norb,nspin))
    allocate(em(norb,2), ek(nbath,2), vmk(norb,nbath,2))

    ! initial impurity hamiltonian parameters
    call initial_h_imp_params( em, ek, vmk )

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

        if (iloop>1) then
            ! find new em, ek, vmk 
            call fit_h_imp_params( em, ek, vmk, D_cl, G_cl, G_loc )
        endif

        call cluster_hybridization_ftn( ek, vmk, D_cl )

        call ground_state( em, ek, vmk, gs )

        call cluster_green_ftn( em, ek, vmk, gs, G_cl )

        call local_green_ftn( em, D_cl, G_cl, G_loc )

        call dump_data( em, ek, vmk, D_cl, G_cl, G_loc )

        call test_convergence( G_cl, G_loc, diff, converged )

        if (master) then
            write(*,*) "DMFT SCF condition |G_cl - G_loc|^2 = ", diff
        endif

        if (converged) then
            if (master) then
                write(*,*) "DMFT loop has converged within ",iloop," iterations."
            endif
            exit dmftloop
        endif

        t2_loop = mpi_wtime(mpierr)
    enddo dmftloop

    call post_process( em, ek, vmk, D_cl, G_cl, G_loc, gs )

    deallocate(G_cl, G_loc, D_cl)
    deallocate(em, ek, vmk)
    deallocate(gs%vec)
end subroutine dmft

subroutine test_convergence( G_cl, G_loc, diff, converged )
    use mpi
    use dmft_params, only: nspin, norb, nbath, nw, scf_tol
    use dmft_grid, only: nwloc

    double complex, intent(in) :: &
        G_cl(nwloc, norb, nspin), &
        G_loc(nwloc, norb, nspin)

    logical, intent(out) :: converged
    double precision :: d, diff_loc, diff

    integer :: iw, iorb, ispin

    diff_loc = 0.d0
    do ispin=1,nspin
        do iorb=1,norb
            do iw=1,nwloc
                d = abs(G_cl(iw,iorb,ispin)-G_loc(iw,iorb,ispin))
                diff_loc = diff_loc + d*d
            enddo
        enddo
    enddo

    diff_loc = diff_loc / (norb*nbath*nw)
    call mpi_allreduce( diff_loc, diff, 1, mpi_double_precision, &
                        mpi_sum, comm, mpierr )

    converged = diff < scf_tol

end subroutine test_convergence
