subroutine g02_rotation_Gamma
!rotation_Amat‚Æ“¯‚¶@‡@‘@ˆÛ•ûŒü‚²‚Æ‚Éalpha‚ğ‰ñ“] ‡A—v‘f‚ÌŒX‚«‚²‚Æ‚Éalpha‚ğ‰ñ“] ŠÈ’P‚É‚µ‚Ü‚µ‚½2023/12/16
use parameter1

implicit none

integer a,b, e, f, i, j, k, l, m, n, p, q,ii,jj
integer i1, i2,  j1, j2, k1,k2

real(8) temp(3,3)!‚È‚ñ‚±‚êHH
allocate(Gamma_mat(3,3))
allocate(rot_Gamma_mat(3,3))
allocate(rot_Gamma_ini(3,3))
allocate(g_rot_Gamma(g_econst,6))
allocate(Qtemp(3,3,4))
allocate(Qrot_Gamma_ini(3,3,4))








do a=1,g_econst
temp=0.0d0
!ƒxƒNƒgƒ‹‚ğ2ŠK‚Ìƒeƒ“ƒ\ƒ‹‚É‘‚«’¼‚µ‚Ä‰ñ“]‚³‚¹‚é
Gamma_mat=0.0d0
do i=1,3
Gamma_mat(i,i)=gamma(i)
end do     !‚¹‚ñ’f¬•ª‚ğ’²®
  Gamma_mat(1,2) = 0.5d0 * gamma(4)
  Gamma_mat(2,1) = 0.5d0 * gamma(4)
  Gamma_mat(2,3) = 0.5d0 * gamma(5)
  Gamma_mat(3,2) = 0.5d0 * gamma(5)
  Gamma_mat(1,3) = 0.5d0 * gamma(6)
  Gamma_mat(3,1) = 0.5d0 * gamma(6)

!‡@Ï‘w•ûŒü‚ğl—¶‚µ‚Ä1²ü‚è‚Ì‰ñ“]iÅ‰‚Í90degj
g_Qmat=0.0d0

g_Qmat(1,1) = 1
g_Qmat(2,2) = cos(deg_f(fiber_dir(a)))
g_Qmat(2,3) = -sin(deg_f(fiber_dir(a)))
g_Qmat(3,2) = sin(deg_f(fiber_dir(a)))
g_Qmat(3,3) = cos(deg_f(fiber_dir(a)))

do i=1,3
do j=1,3
do k=1,3
temp(i,j)=temp(i,j)+g_Qmat(k,i)*Gamma_mat(k,j)   
end do
end do
end do


rot_Gamma_ini=0.0d0
do i=1,3
do j=1,3
do k=1,3
rot_Gamma_ini(i,j)=rot_Gamma_ini(i,j)+temp(i,k)*g_Qmat(k,j)
end do
end do
end do



!‡A—v‘f‚ÌŒX‚«‚É‰‚¶‚½3²ü‚è‚Ì‰ñ“]
  g_Qmat_2=0.0d0
 g_Qmat_2(1,1)=cos(-th(a))
 g_Qmat_2(1,2)=-sin(-th(a))
 g_Qmat_2(2,1)=sin(-th(a))
 g_Qmat_2(2,2)=cos(-th(a))
 g_Qmat_2(3,3)=1

Qtemp=0.0d0

do i=1,3
do j=1,3
do k=1,3
Qtemp(i,j,fiber_dir(a))=Qtemp(i,j,fiber_dir(a))+g_Qmat_2(k,i)*rot_Gamma_ini(k,j)
end do
end do
end do
Qrot_Gamma_ini=0.0d0
do i=1,3
do j=1,3
do k=1,3
Qrot_Gamma_ini(i,j,fiber_dir(a))=Qrot_Gamma_ini(i,j,fiber_dir(a))+Qtemp(i,k,fiber_dir(a))*g_Qmat_2(k,j)  !‚±‚Ì‚ ‚½‚èfiber_dir‚Ì€‚Í—v‚ç‚È‚¢‹C‚ª‚·‚é@‚»‚Ì‚¤‚¿‚â‚é@2023/12/16
end do
end do
end do


!Ši”[
do i=1,3
   g_rot_Gamma(a,i) = Qrot_Gamma_ini(i,i,fiber_dir(a))
  end do
  g_rot_Gamma(a,4) = 2.0d0 * Qrot_Gamma_ini(1,2,fiber_dir(a)) 
  g_rot_Gamma(a,5) = 2.0d0 * Qrot_Gamma_ini(2,3,fiber_dir(a))
  g_rot_Gamma(a,6) = 2.0d0 * Qrot_Gamma_ini(3,1,fiber_dir(a))

end do


open(10,file='Output_rot_Gamma.dat')
  do a=1,g_econst
write(10,'("element",I5)')a
      write(10,'(6F15.10)') (g_rot_Gamma(a, i), i=1,6)
    end do

close(10)

write(*,*) "rotation_Gamma COMPLETE!!"




end subroutine