module ed_params
!===============================================================================
! This module contains the option parameters used only in ed_solver module.
! Note that parameters declared in this module must not be changed after
! the initialization. 
!===============================================================================

    use dmft_params
    use mpi
    use numeric_utils, only: icom

    implicit none

    public :: ed_read_params

    integer, public :: &
        nbathperorb, &   ! nbath per orbital (when orbital diagonal)

    double precision, allocatable, public :: &
        ek_in(:),    &   ! ek_in(nsite)       initial impurity/bath levels
        vk_in(:,:)       ! vk_in(norb,nbath)  initial impurity-bath hybridization

    character(len=100), public :: diag_method

    logical, public :: &
        print_arpack_stat,  &
        read_h_imp_params,  & ! read initial impurity hamiltonian parameters 
                              ! from file h_imp.params
        em_present, ek_present, vk_present

    private
contains

    subroutine ed_read_params
        use fdf
        use utils
        character(len=200) :: text
        character(len=100), parameter :: &
            FMT_INT    = "(3x,a40,2x,a,2x,I8)",   &
            FMT_DOUBLE = "(3x,a40,2x,a,2x,F8.3)", &
            FMT_EXP    = "(3x,a40,2x,a,2x,ES8.1)"

        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline
        integer :: i,j,ispin

        diag_method = fdf_get("DMFT.ED.Diagonalization", "full") 

        if (mod(nbath,norb)/=0) then
            call die("ed_read_params", "nbath should be multiple of norb.")
        endif
        nbathperorb = nbath/norb
        nsite = norb+nbath
        select case(diag_method)
            case ("full")
                if (nsite>8) then
                    call die("ed_read_params", "full diagonalization for nsite > 8 is not recommended. use arpack.")
                endif
            case default

        end select

        read_h_imp_params = fdf_get("DMFT.ED.HamiltonianParamsFromFile", .false.)

        allocate(ek_in(nsite),vk_in(norb,nbath))
        ek_in = 0.0d0
        vk_in = 0.0d0
        em_present = fdf_block('DMFT.ED.InitialImpurityLevels', bfdf)
        if (em_present) then
            i = 1
            do while( (i .le. norb) .and. (fdf_bline(bfdf, pline)))
                ek_in(i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        ek_present = fdf_block('DMFT.ED.InitialBathLevels', bfdf)
        if (fdf_block('DMFT.ED.InitialBathLevels', bfdf)) then
            i = 1
            do while( (i .le. nbath) .and. (fdf_bline(bfdf, pline)))
                ek_in(Norb+i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        vk_present = fdf_block('DMFT.ED.InitialHybridization', bfdf)
        if (vk_present) then
            i = 1
            do while( i.le.nbath .and. (fdf_bline(bfdf, pline)))
                do j=1,norb
                    vk_in(j,i) = fdf_breals(pline,j)
                enddo
                i = i + 1
            enddo
        endif


        print_arpack_stat = fdf_get("DMFT.ED.ARPACKStat", .false.)

        if (master) then
            write(6,*)
            write(6,'(a)') repeat("=",80)
            write(6,'(4x,a)') "Exact Diagoanlization Solver Parameters"
            write(6,'(a)') repeat("=",80)
            write(text,*) "Diagonalization Method"
            write(6,'(3x,a40,2x,a,2x,a)') text, '=', diag_method
            write(text,*) "Number of impurity orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Norb
            write(text,*) "Number of bath orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nbath
            write(text,*) "Number of sites"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nsite
            write(text,*) "Number of eigenvalues to be computed"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', nev
            write(text,*) "Number of steps in continued fraction"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nstep
            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif
    end subroutine ed_read_params
end module ed_params
