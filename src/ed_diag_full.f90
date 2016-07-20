module ed_diag_full
!==============================================================================
! Serial, full diagonalization for each sector using LAPACK.
! This is for testing small problems.
! Currently, only one sector can be diagonalized in a run.
!==============================================================================
    
    use mpi
    use dmft_params, only: norb, beta
    use ed_params, only: nsite, nbath, sectors, nsector, PROB_THRESHOLD, eigpair_t
    use ed_basis, only: generate_basis, basis_t, ed_basis_get
    use ed_hamiltonian, only: generate_hamiltonian, ek,vk
    use numeric_utils, only: boltzmann_factor
    use utils, only: die

    implicit none

    public :: diag_full

    private
contains

    subroutine diag_full(ia, nev_calc, eigpairs)
        integer, intent(in) :: ia
        integer, intent(out) :: nev_calc
        type(eigpair_t), allocatable, intent(out) :: eigpairs(:)

        ! local variables
        integer :: isector, iev, ne_up, ne_down, nloc, nh
        double precision, allocatable :: H(:,:), ev(:), prob(:)
        double precision :: Z
        type(basis_t) :: basis

        integer :: i,j

        if (nsector.gt.1) then
            call die("diag_full", &
                "full diagonalization for nsector>1 is not implemented.")
        else if (nprocs.gt.1) then
            call die("diag_full", &
                "full diagonalization for nprocs>1 is not implemented.")
        endif

        if (master) then
            write(*,"(x,a,I3)") "Full diagonalization..."
        endif

        ne_up = sectors(1,1)
        ne_down = sectors(1,2)
        nh = sectors(1,3)

        ! basis states for the sector (ne_up,ne_down)
        call generate_basis(ne_up, ne_down, basis)
            
        allocate(H(nh,nh),ev(nh),prob(nh))

        call generate_hamiltonian(ia,basis,H)

        ! note that lapack returned eigenvalues will be in ascending order
        call lapack_diag(basis,H,ev)

        nev_calc = 1
        prob(1) = 1.0d0 ! ground state
        do iev=2,nh
            prob(iev) = boltzmann_factor(beta, ev(iev)-ev(1))
            if (prob(iev).gt.PROB_THRESHOLD) then
                nev_calc = nev_calc + 1
            endif
        enddo
        Z = sum(prob(1:nev_calc))

        allocate(eigpairs(nev_calc))

        do iev=1,nev_calc
            eigpairs(iev)%sector = 1
            eigpairs(iev)%level  = iev
            eigpairs(iev)%val    = ev(iev)
            eigpairs(iev)%prob   = prob(iev)/Z
            allocate(eigpairs(iev)%vec(sectors(1,3)))
            eigpairs(iev)%vec(:) = H(:,iev)
        enddo

        if (master) then
            write(6,*)
            write(6,*) "Obatined eigenvalues."
            write(6,"(a,ES10.3,a,I5)") " Number of eigenvalues (with prob > ", &
                                       PROB_THRESHOLD,") = ", nev_calc
            write(6,*)
            write(6,"(a)") " Eigenvalue          Prob           Sector   Level"
            do iev=1,nev_calc
                write(6,"(1x,F16.10,4x,ES13.5,2I8)") eigpairs(iev)%val,eigpairs(iev)%prob,&
                                          eigpairs(iev)%sector,eigpairs(iev)%level
            enddo
            write(6,*)
        endif
    end subroutine diag_full

    subroutine lapack_diag(basis,H,ev)
        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: H(basis%ntot,basis%ntot)
        double precision, intent(out) :: ev(basis%ntot)
        ! local variables
        integer :: nh, info
        double precision :: work(3*basis%ntot)

        nh = basis%ntot

        call DSYEV('V','U',basis%ntot,H,basis%ntot,ev,work,3*basis%ntot,info)

        if (info.ne.0) then
            write(*,*) "DSYEV info = ", info
            call die("lapack_diag","failed to diagonalize")
        endif

    end subroutine lapack_diag

end module ed_diag_full
