module ed_eigpair

    implicit none

    type, public :: eigpair_t
        integer :: sector  ! sector index
        integer :: level   ! level index within the sector
        double precision :: val   ! eigenvalue 
        double precision :: prob  ! exp(-beta*val)/Z

        integer :: nloc    ! dimension of the eigenvector local to the processor
        double precision, pointer :: vec(:) ! eigenvector
    end type eigpair_t

    ! number of eigenvalues that has prob>PROB_THRESHOLD
    ! and corresponding eigenvalues & eigenvectors
    integer :: nev_calc
    type(eigpair_t), allocatable :: eigpairs(:)

contains

    subroutine allocate_eigpairs(nev)
        integer, intent(in) :: nev

        allocate(eigpairs(nev))

    end subroutine allocate_eigpairs

end module ed_eigpair
