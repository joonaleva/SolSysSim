module vectors
    use constants
    implicit none
    type :: vector                  ! define the vector type, an object containing 3 real numbers: x, y and z
        real(rk) :: x,y,z
    end type vector

    interface operator(+)           ! define (+) operator as a vector addition
        module procedure vecadd             ! vec + vec
    end interface operator(+)

    interface operator(-)           ! define (-) operator as a vector subtraction
        module procedure vecsub             ! vec - vec
        module procedure negvec             ! -vec
    end interface operator(-)

    interface operator(*)           ! define (*) operator as a product of a constant and a vector, and the dot product of two vectors
        module procedure vecprod_left       ! const * vec
        module procedure vecprod_right      ! vec * const
        module procedure dotp               ! vec * vec 
    end interface operator(*)

    interface operator(/)           ! define (//) operator as a division of a vector and a constant
        module procedure vecdiv             ! vec / const
    end interface operator(/)

    contains
    elemental function vecadd(vec1, vec2) result(res) !addition of vectors
        implicit none
        type(vector),intent(in) :: vec1, vec2
        type(vector) :: res

        res%x = vec1%x + vec2%x
        res%y = vec1%y + vec2%y
        res%z = vec1%z + vec2%z
    end function 

    function vecsum(vecarr) result(res)   !sum of vector array
        implicit none
        type(vector),intent(in) :: vecarr(:)
        type(vector) :: res
        integer :: i

        res = vector(0, 0, 0)
        do i=1, size(vecarr)
            res = res + vecarr(i)
        end do
    end function

    elemental function vecsub(vec1, vec2) result(res) !subtraction of vectors
        implicit none
        type(vector),intent(in) :: vec1, vec2
        type(vector) :: res

        res%x = vec1%x - vec2%x
        res%y = vec1%y - vec2%y
        res%z = vec1%z - vec2%z
    end function

    elemental function negvec(vec) result(res) !negative of a vector
        implicit none
        type(vector),intent(in) :: vec
        type(vector) :: res

        res%x = -vec%x
        res%y = -vec%y
        res%z = -vec%z
    end function

    elemental function vecprod_left(const,vec) result(res) !product constant and vector
        implicit none
        real(rk),intent(in) :: const
        type(vector),intent(in) :: vec
        type(vector) :: res

        res%x = const * vec%x
        res%y = const * vec%y
        res%z = const * vec%z
    end function

    elemental function vecprod_right(vec,const) result(res) !product of constant and vector
        implicit none
        type(vector),intent(in) :: vec
        real(rk),intent(in) :: const
        type(vector) :: res

        res%x = vec%x * const
        res%y = vec%y * const
        res%z = vec%z * const
    end function

    elemental function vecdiv(vec, const) result(res) !division of a vector with a constant
        implicit none
        real(rk),intent(in) :: const
        type(vector),intent(in) :: vec
        type(vector) :: res

        res%x = vec%x / const
        res%y = vec%y / const
        res%z = vec%z / const
    end function

    elemental function dotp(vec1, vec2) result(res) !dot product of vectors
        implicit none
        type(vector),intent(in) :: vec1, vec2
        real(rk) :: res

        res = vec1%x * vec2%x + vec1%y * vec2%y + vec1%z * vec2%z
    end function 

    elemental function crossp(vec1, vec2) result(res) !cross product of vectors
        implicit none
        type(vector),intent(in) :: vec1, vec2
        type(vector) :: res

        res%x = vec1%y*vec2%z - vec1%z*vec2%y
        res%y = vec1%z*vec2%x - vec1%x*vec2%z
        res%z = vec1%x*vec2%y - vec1%y*vec2%x
    end function

    elemental function veclen(vec) result(res) !length of a vector
        implicit none
        type(vector),intent(in) :: vec
        real(rk) :: res

        res = sqrt(vec%x**2 + vec%y**2 + vec%z**2)
    end function

    elemental function vecdist(vec1, vec2) result(res) !distance between 2 points given with vectors
        implicit none
        type(vector),intent(in) :: vec1, vec2
        real(rk) :: res

        res = veclen(vec1 - vec2)
    end function

    elemental function unitvec(vec) result(res) !unit vector of given vector
        implicit none
        type(vector),intent(in) :: vec
        type(vector) :: res
        real(rk) :: L

        L = veclen(vec)
        if (L==0.0_rk) then             ! failsafe for vectors (0, 0, 0) so no division by 0 occurs
            res = vector(0, 0, 0)
        else
            res = vec / L
        end if
    end function

end module vectors