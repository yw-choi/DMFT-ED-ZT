module io_units

    implicit none

    integer, parameter :: IO_TB_DOS = 11
    character(len=100), parameter :: FN_TB_DOS = "tbdos.dat"

    integer, parameter :: IO_GR     = 100
    character(len=100), parameter :: FN_GR = "green.dat"

    integer, parameter :: IO_H_PARAMS = 101
    character(len=100), parameter :: FN_H_PARAMS = "h_imp.dat"

    integer, parameter :: IO_SIGMA = 102
    character(len=100), parameter :: FN_SIGMA = "sigma.dat"

    integer, parameter :: IO_G_COEFFS = 103
    character(len=100), parameter :: FN_G_COEFFS = "g_coeffs.dat"
end module io_units
