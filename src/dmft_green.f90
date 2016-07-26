module dmft_green
    use mpi
    use io_units
    use dmft_params, only: norb, nbath, nspin, sectors, maxnstep
    use dmft_grid, only: nwloc, omega
    use ed_basis, only: basis_t, dealloc_basis, generate_basis
    use eigpair, only: eigpair_t
    use constants

    implicit none

    public :: &
        cluster_hybridization_ftn, &
        cluster_green_ftn, &
        local_green_ftn

    private 
contains
    ! =========================================================================
    ! PUBLIC subroutines
    ! =========================================================================

    ! calculates the cluster hybridization function from ek, vmk 
    subroutine cluster_hybridization_ftn( ek, vmk, D_cl )
        double precision, intent(in) :: &
            ek(nbath,2), vmk(norb,nbath,2)
        double complex, intent(out) :: &
            D_cl(nwloc,norb,nspin)

        integer :: iw,iorb,ispin,ibath
        double complex :: d, iwn

        do ispin=1,nspin
            do iorb=1,norb
                do iw=1,nwloc
                    D_cl(iw,iorb,ispin) = 0.d0
                    do ibath=1,nbath
                        iwn = cmplx(0.d0,omega(iw)) 
                        d = vmk(iorb,ibath,ispin)**2/(iwn-ek(ibath,ispin))
                        D_cl(iw,iorb,ispin) = D_cl(iw,iorb,ispin) + d
                    enddo
                enddo
            enddo
        enddo
    end subroutine cluster_hybridization_ftn 

    subroutine cluster_green_ftn( em, ek, vmk, gs, G_cl )
        double precision, intent(in) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        type(eigpair_t), intent(in) :: gs
        
        double complex, intent(out) :: G_cl(nwloc,norb,nspin)

        logical :: even
        integer :: isector, ne_u, ne_d, iorb
        type(basis_t) :: basis
        double complex, allocatable :: g_diag(:,:)

        G_cl = 0.d0
        allocate(g_diag(nwloc,2))

        do iorb=1,norb
            isector = gs%isector
            ne_u = sectors(isector,1)
            ne_d = sectors(isector,2)

            if (nspin==2.or.ne_u/=ne_d) then
                even = .false.
            else
                even = .true.
            endif

            call generate_basis(ne_u, ne_d, basis)
            
            ! G_{iorb,up}
            call green_diag( iorb, 1, basis, gs, em, ek, vmk, &
                             g_diag )

            if (.not.even) then
                ! if not spin even, calculates G_{iorb,dn}
                call green_diag( iorb, 2, basis, gs, em, ek, vmk, &
                                 g_diag )
            endif

            if (nspin==1 .and. even) then
                G_cl(:,iorb,1) = g_diag(:,1)
            else if (nspin==1 .and. .not.even) then
                G_cl(:,iorb,1) = 0.5d0*(g_diag(:,1)+g_diag(:,2))
            else ! nspin==2
                G_cl(:,iorb,1) = g_diag(:,1)
                G_cl(:,iorb,2) = g_diag(:,2)
            endif

            call dealloc_basis(basis)
        enddo

    end subroutine cluster_green_ftn 

    subroutine local_green_ftn( em, D_cl, G_cl, G_loc )
        double precision, intent(in) :: &
            em(norb,nspin)

        double complex, intent(in) :: &
            G_cl(nwloc,norb,nspin), &
            D_cl(nwloc,norb,nspin)

        double complex, intent(out) :: G_loc(nwloc,norb,nspin)

        integer :: iw, iorb, ispin, ix, ne

        double precision :: val,x

        ne = 10000

        ! @TODO implemented only for the Bethe-lattice
        ! G_loc = 4.d0*D_cl
        do ispin=1,nspin
            do iorb=1,norb
                do iw=1,nwloc
                    G_loc(iw,iorb,ispin) = 0.d0

                    do ix=1,ne
                        x = -1.d0+(ix-1)*2.d0/(ne-1)
                        val = 1.d0-x*x
                        G_loc(iw,iorb,ispin) = G_loc(iw,iorb,ispin) + &
                            2.d0/pi*sqrt(max(val,0.d0))/&
                            (em(iorb,ispin) &
                            + D_cl(iw,iorb,ispin) &
                            + 1.0d0/G_cl(iw,iorb,ispin) &
                            - x)
                    enddo
                    G_loc(iw,iorb,ispin) = G_loc(iw,iorb,ispin)/(ne-1)*2.d0
                enddo
            enddo
        enddo
        
    end subroutine local_green_ftn


    ! =========================================================================
    ! PRIVATE subroutines
    ! =========================================================================

    subroutine green_diag( iorb, ispin, basis, gs, em, ek, vmk, g_diag )
        ! @TODO move lanczos_iteration subroutine to a separate module
        use lanczos_solver, only: lanczos_iteration
        use numeric_utils, only: continued_fraction_m, continued_fraction_p, &
                                 mpi_norm
        use ed_operator, only: apply_c

        integer, intent(in) :: iorb, ispin
        type(basis_t), intent(in) :: basis
        type(eigpair_t), intent(in) :: gs
        double precision, intent(in) :: em(norb,2), ek(nbath,2), vmk(norb,nbath,2)
        double complex, intent(out) :: g_diag(nwloc,2)

        type(basis_t) :: basis_out
        integer :: nstep, nlocup, iw
        double precision, allocatable :: v_init(:), a(:), b(:)
        double complex :: z, gr

        allocate(a(maxnstep),b(maxnstep))

        g_diag(:,ispin) = 0.d0

        ! G^+
        ! 1. c^+|iev>
        if (ispin == 1) then
            call generate_basis( basis%ne_up+1, basis%ne_down, &
                                 basis_out)
        else
            call generate_basis( basis%ne_up, basis%ne_down+1, &
                                 basis_out)
        endif

        allocate(v_init(basis_out%nloc))
        nlocup = basis_out%nloc
        
        call apply_c(basis, gs%vec, 1, iorb, ispin, basis_out, v_init)

        ! 2. lanczos coefficients
        call lanczos_iteration( em, ek, vmk, basis_out, v_init, nstep, a, b )
        b(1) = mpi_norm(v_init, basis_out%nloc)

        ! 3. green function as a continued fraction
        do iw=1,nwloc
            z = cmplx(gs%val, omega(iw))
            gr = continued_fraction_p(z, nstep, a, b)
            g_diag(iw,ispin) = g_diag(iw,ispin)+gr
        enddo

        call dealloc_basis(basis_out)

        ! G^-
        ! 1. c^-_{iorb,ispin} | eigvec >
        if (ispin == 1) then
            call generate_basis( basis%ne_up-1, basis%ne_down, &
                                 basis_out)
        else
            call generate_basis( basis%ne_up, basis%ne_down-1, &
                                 basis_out)
        endif

        if (nlocup /= basis_out%nloc) then
            deallocate(v_init)
            allocate(v_init(basis_out%nloc))
        endif

        call apply_c(basis, gs%vec, 2, iorb, ispin, basis_out, v_init)

        ! 2. lanczos coefficients
        call lanczos_iteration( em, ek, vmk, basis_out, v_init, nstep, a, b )
        b(1) = mpi_norm(v_init, basis_out%nloc)

        ! 3. green function as a continued fraction
        do iw = 1,nwloc
            z = cmplx(-gs%val, omega(iw))
            gr = continued_fraction_m(z, nstep, a, b)
            
            g_diag(iw,ispin) = g_diag(iw,ispin)+gr
        enddo

        call dealloc_basis(basis_out)
    end subroutine green_diag
end module dmft_green
