module funcs_f
    contains
        real*8 function f1 (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
              & 3.0 * e**2 + &
              & 0.375 * e**4) / &
              & (1. - e**2)**(4.5)
        end function f1

        real*8 function f2 (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
              & 7.5 * e**2 + &
              & 5.625 * e**4 + &
              & 0.3125 * e**6) / &
              & (1. - e**2)**6
        end function f2
        
        real*8 function f3 (e) result (f)
        implicit none
        real*8, intent (in) :: e
            f = (1. + &
              & 15.5 * e**2 + &
              & 31.875 * e**4 + &
              & 11.5625 * e**6 + &
              & 0.390625 * e**8) / &
              & (1. - e**2)**(7.5)
        end function f3
        
        real*8 function f4 (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
              & 1.5 * e**2 + &
              & 0.125 * e**4) / &
              & (1. - e**2)**5
        end function f4
        
        real*8 function f5 (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = (1. + &
              & 3.75 * e**2 + &
              & 1.875 * e**4 + &
              & 0.078125 * e**6) / &
              & (1. - e**2)**(6.5)
        end function f5

        real*8 function f1_short (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 7.5 * e**2
        end function f1_short

        real*8 function f2_short (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 13.5 * e**2
        end function f2_short

        real*8 function f3_short (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 23. * e**2
        end function f3_short

        real*8 function f4_short (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 6.5 * e**2
        end function f4_short

        real*8 function f5_short (e) result (f)
            implicit none
            real*8, intent (in) :: e
            f = 1. + 10.25 * e**2
        end function f5_short
end module funcs_f