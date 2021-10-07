module f_funcs
    
    implicit none

    abstract interface 

        real*8 function f_tem (e) result (f)
            real*8 :: e
        end function f_tem

    end interface

    contains

        real*8 function f1_l (e) result (f)
            implicit none
            real*8, intent (in) :: e
                f = (1. + &
                & 3.0 * e**2 + &
                & 0.375 * e**4) / &
                & (1. - e**2)**(4.5)
        end function f1_l

        real*8 function f2_l (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
            & 7.5 * e**2 + &
            & 5.625 * e**4 + &
            & 0.3125 * e**6) / &
            & (1. - e**2)**6
        end function f2_l
        
        real*8 function f3_l (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
            & 15.5 * e**2 + &
            & 31.875 * e**4 + &
            & 11.5625 * e**6 + &
            & 0.390625 * e**8) / &
            & (1. - e**2)**(7.5)
        end function f3_l
        
        real*8 function f4_l (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
            & 1.5 * e**2 + &
            & 0.125 * e**4) / &
            & (1. - e**2)**5
        end function f4_l
        
        real*8 function f5_l (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
            & 3.75 * e**2 + &
            & 1.875 * e**4 + &
            & 0.078125 * e**6) / &
            & (1. - e**2)**(6.5)
        end function f5_l

        real*8 function f1_s (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 7.5 * e**2
        end function f1_s

        real*8 function f2_s (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 13.5 * e**2
        end function f2_s

        real*8 function f3_s (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 23. * e**2
        end function f3_s

        real*8 function f4_s (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 6.5 * e**2
        end function f4_s

        real*8 function f5_s (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 10.25 * e**2
        end function f5_s

end module f_funcs