subroutine solve_thermal_elastic

use parameter1

implicit none

integer i,j,k,l,m,n

call make_BTgamma
call solve_disp
call solve_strain
call solve_stress
!Ç±ÇÍÇ≈èIÇÌÇËÇ≈Ç¶Ç¶ÇÒÇ©Ç‚ÅH

end subroutine