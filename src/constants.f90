module constants
    implicit none
    integer, parameter :: ik = 4, rk = 8    ! integer and real kinds

    !integer(), parameter :: au = 149597870700  !meters, integer size exceeded

    real(rk),parameter :: pi = 4*atan(1.0_rk)
    real(rk),parameter :: G_SI = 6.67430e-11_rk !Newtonian gravitation constant in SI units, m^3 kg^-1 s^-2
    real(rk),parameter :: day = 86400           !seconds
    real(rk),parameter :: year = 365.25_rk*day  !Julian year in seconds
    real(rk),parameter :: au = 1.495978707e11_rk    !astronomical unit in meters
    real(rk),parameter :: pc = 6.48e5_rk/pi*au      !parsec in meters
    real(rk),parameter :: M_sol = 1.988416e30_rk    !Solar mass in kilograms

    !real(rk),parameter :: G_sol = 4.3009172706e-3_rk  !pc M_sol^-1 (km/s)^2
    real(rk),parameter :: G_aud = G_SI / (au**3*day**(-2))  !Newtonian gravitation constant in au^3 kg^-1 day^-2
    real(rk),parameter :: G_sol = G_SI / (pc**3 * M_sol**(-1) * year**(-2))    !Newtonian gravitation constant in pc^3 M_sol^-1 year^
end module constants