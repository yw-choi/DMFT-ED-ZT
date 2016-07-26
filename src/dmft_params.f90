module dmft_params
! =============================================================================
! DMFT parameters that do not change during the run.
! =============================================================================
    implicit none

    integer ::          &
        nspin,          & ! number of spin components 
        norb,           & ! number of impurity orbitals
        nbath,          & ! number of bath sites
        nsite,          & ! norb+nbath
        nloop,          & ! maximum number of DMFT loop
        nsector,        & ! number of hamiltonain sector in (Q,Sz) basis
        maxnstep,       & ! number of steps in calculating continued fraction 
        nw,             & ! Total number of matsubara frequencies
        diag_solver,    & ! diagonalization solver
        tbham             ! type of tight-binding hamiltonian
                          ! 0 : read from file
                          ! 1 : 2d square lattice
                          ! 2 : bethe lattice, infinite dimension (circular dos)

    logical :: &
        em_present,         &
        ek_present,         &
        vmk_present,        &
        read_gf_from_file

    double precision ::     &
        U,                  & ! U
        Up,                 & ! U', U'=U-2J with rotational invariance
        Jex,                & ! J 
        Jp,                 & ! J' 
        Mu,                 & ! chemical potential
        beta,               & ! fictitious inverse temperature (low-frequency cutoff)
        scf_tol               ! DMFT SCF tolerance

    double precision, allocatable :: &
        em_input(:,:),       & ! em_input(norb,2)        initial impurity levels
        ek_input(:,:),       & ! em_input(norb,2)        initial bath levels
        vmk_input(:,:,:)       ! vmk_input(norb,nbath,2) initial hybridizations

    integer, allocatable, public :: &
        sectors(:,:)  ! sectors(nsector,3) ne_up/ne_down/nbasis in each sectors

contains

    ! Read general DMFT parameters independent of solver.
    subroutine dmft_params_read
        use mpi
        use fdf
        use numeric_utils, only: icom
        use utils, only: die

        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline
        character(len=200) :: msg
        character(len=100), parameter ::          &
            FMT_INT     = "(3x,a40,2x,a,2x,I8)",   &
            FMT_LOGICAL = "(3x,a40,2x,a,2x,L8)",   &
            FMT_DOUBLE  = "(3x,a40,2x,a,2x,F8.3)", &
            FMT_EXP     = "(3x,a40,2x,a,2x,ES8.1)"
        integer :: i,j,k

        if (master) then
            write(6,*)
            write(6,'(a)') repeat("=",80)
            write(6,'(4x,a)') "DMFT Parameters"
            write(6,'(a)') repeat("=",80)
        endif

        nspin = fdf_get("Nspin",1)
        norb  = fdf_get("Norb",1)
        nbath = fdf_get("Nbath",5)
        nsite = norb + nbath

        U   = fdf_get("U",   1.0D0)
        Jex = fdf_get("J",   0.3D0) 
        Up  = fdf_get("Up", -1.0d0)
        Jp  = fdf_get("Jp", -1.0d0) 

        ! rotational invariance
        if (Up<0) then
            Up = U - 2*Jex
        endif

        if (Jp<0) then
            Jp = Jex
        endif

        Mu  = fdf_get("Mu",0.5d0)
        beta = fdf_get("beta", 200)

        nloop = fdf_get("MaxDMFTIterations", 100)
        scf_tol = fdf_get("SCFTolerance", 1.d-6)
        maxnstep = fdf_integer("LanczosSteps", 80)

        nw = fdf_get("Nw", 1000)

        tbham = fdf_get("TBHamiltonian", 2)

        nsector = fdf_get("nsector", 1)
        allocate(sectors(nsector,3))

        if (fdf_block('sectors', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nsector))
                sectors(i,1) = fdf_bintegers(pline,1) ! nup
                sectors(i,2) = fdf_bintegers(pline,2) ! ndown

                ! sector dimension
                sectors(i,3) = icom(sectors(i,1)+sectors(i,2),sectors(i,1))* &
                               icom(sectors(i,1)+sectors(i,2),sectors(i,2))
                            
                i = i + 1
            enddo
        endif

        allocate(em_input(norb,2), ek_input(nbath,2))
        allocate(vmk_input(norb,nbath,2))

        em_present = fdf_block('InitialImpurityLevels', bfdf)
        if (em_present) then
            do i=1,norb
                if (fdf_bline(bfdf, pline)) then
                    em_input(i,:) = fdf_breals(pline, 1)
                    if (nspin==2) then
                        em_input(i,2) = fdf_breals(pline, 2)
                    endif
                else
                    call die("dmft_params", &
                        "not enough initial impurity levels given")
                endif
            enddo
        endif

        ek_present = fdf_block('InitialBathLevels', bfdf)
        if (ek_present) then
            do i=1,nbath
                if (fdf_bline(bfdf, pline)) then
                    ek_input(i,:) = fdf_breals(pline, 1)
                    if (nspin==2) then
                        ek_input(i,2) = fdf_breals(pline, 2)
                    endif
                else
                    call die("dmft_params", &
                        "not enough initial bath levels given")
                endif
            enddo
        endif

        vmk_present = fdf_block('InitialHybridization', bfdf)
        if (vmk_present) then
            do i=1,nbath
                if (fdf_bline(bfdf, pline)) then
                    do j=1,norb
                        vmk_input(j,i,:) = fdf_breals(pline, j)
                    enddo
                else
                    call die("dmft_params", &
                        "not enough initial hybridization parameters given")
                endif
            enddo

            if (nspin==2) then
                do i=1,nbath
                    if (fdf_bline(bfdf, pline)) then
                        do j=1,norb
                            vmk_input(j,i,2) = fdf_breals(pline, j)
                        enddo
                    else
                        call die("dmft_params", &
                            "not enough initial hybridization parameters given")
                    endif
                enddo
            endif
        endif

        read_gf_from_file = fdf_get("ReadGreenFtn", .false.)

        diag_solver = fdf_get("DiagSolver", 1)

        if (master) then
            write(msg,*) "Number of impurity orbitals"
            write(6,FMT_INT) msg, '=', Norb

            write(msg,*) "Number of bath orbitals"
            write(6,FMT_INT) msg, '=', Nbath

            write(msg,*) "Number of spin components"
            write(6,FMT_INT) msg, '=', nspin

            write(msg,*) "Number of sites"
            write(6,FMT_INT) msg, '=', nsite

            write(msg,*) "U"
            write(6,FMT_DOUBLE) msg, '=', U

            write(msg,*) "Up"
            write(6,FMT_DOUBLE) msg, '=', Up

            write(msg,*) "J"
            write(6,FMT_DOUBLE) msg, '=', Jex

            write(msg,*) "Jp"
            write(6,FMT_DOUBLE) msg, '=', Jp

            write(msg,*) "chemical potential"
            write(6,FMT_DOUBLE) msg, '=', Mu

            write(msg,*) "Maximum DMFT iterations"
            write(6,FMT_INT) msg, '=', Nloop

            write(msg,*) "DMFT SCF tolerance"
            write(6,FMT_EXP) msg, '=', scf_tol

            write(msg,*) "Number of steps in continued fraction"
            write(6,FMT_INT) msg, '=', maxnstep

            write(msg,*) "Number of Matsubara Frequencies"
            write(6,FMT_INT) msg, '=', nw

            write(msg,*) "Tight-binding Hamiltonian"
            write(6,FMT_INT) msg, '=', tbham

            write(msg,*) "Diagonalization Solver"
            write(6,FMT_INT) msg, '=', diag_solver

            write(msg,*) "Read Green's function from file"
            write(6,FMT_LOGICAL) msg, '=', read_gf_from_file

            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif

    end subroutine dmft_params_read

end module dmft_params
