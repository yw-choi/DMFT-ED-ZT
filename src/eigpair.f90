module eigpair

    implicit none

    type eigpair_t
        integer :: isector
        double precision :: val
        double precision, allocatable :: vec(:)
    end type eigpair_t

end module eigpair
