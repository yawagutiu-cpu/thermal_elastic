subroutine rotation_g_E
!ミクロでの状態を知るためにユニットセルに付加するマクロひずみ(g_E)を正しい方向に回転させる　粘塑性でg_E_speedを回転させるのと一緒
!逆向きに回転させるんじゃよ（　＾ω＾）①要素の傾きを考慮　②積層方向を考慮　の順番？
use parameter1
implicit none

integer a,b,i,j,k,l,p,q
integer i1,i2,j1,j2,k1,k2
real(8) temp(3,3)                             !計算結果を一時的に保存しておくため
allocate(g_E_mat(3,3))
allocate(rot_g_E_ini(3,3))
allocate (temp2(3,3,4))
allocate  (Qrot_g_E_ini(3, 3, 4))




do a=1,g_econst


!①要素の傾きに応じた3軸周りの回転  逆方向回転なので-thでなくth
 g_Qmat_2=0.0d0
 g_Qmat_2(1,1)=cos(th(a))
 g_Qmat_2(1,2)=-sin(th(a))
 g_Qmat_2(2,1)=sin(th(a))
 g_Qmat_2(2,2)=cos(th(a))
 g_Qmat_2(3,3)=1
 
do b=1,g_ipn
g_E_mat=0.0d0
!2階のテンソルに書き直す
do i=1,3
g_E_mat(i,i)=g_E(a,b,i)
end do
g_E_mat(1,2)=0.5d0*g_E(a,b,4)
g_E_mat(2,1)=0.5d0*g_E(a,b,4)
g_E_mat(2,3)=0.5d0*g_E(a,b,5)
g_E_mat(3,2)=0.5d0*g_E(a,b,5)
g_E_mat(1,3)=0.5d0*g_E(a,b,6)
g_E_mat(3,1)=0.5d0*g_E(a,b,6)


temp=0.0d0
do i=1,3
do j=1,3
do k=1,3
temp(i,j)=temp(i,j)+g_Qmat_2(k,i)*g_E_mat(k,j)
end do
end do
end do
rot_g_E_ini=0.0d0
do i=1,3
do j=1,3
do k=1,3
rot_g_E_ini(i,j)=rot_g_E_ini(i,j)+temp(i,k)*g_Qmat_2(k,j)
end do
end do
end do

!②積層方向を考慮して1軸周りの回転（最初は90deg）逆方向回転なのに注意
g_Qmat=0.0d0

g_Qmat(1,1)=1
g_Qmat(2,2)=cos(-deg_f(fiber_dir(a)))
g_Qmat(2,3)=-sin(-deg_f(fiber_dir(a)))
g_Qmat(3,2)=sin(-deg_f(fiber_dir(a)))
g_Qmat(3,3)=cos(-deg_f(fiber_dir(a)))



temp2=0.0d0

do i=1,3
do j=1,3
do k=1,3
temp2(i,j,fiber_dir(a))=temp2(i,j,fiber_dir(a))+g_Qmat(k,i)*rot_g_E_ini(k,j)
end do
end do
end do
Qrot_g_E_ini=0.0d0
do i=1,3
do j=1,3
do k=1,3
Qrot_g_E_ini(i,j,fiber_dir(a))=Qrot_g_E_ini(i,j,fiber_dir(a))+temp2(i,k,fiber_dir(a))*g_Qmat(k,j)
end do
end do
end do



!格納
do i=1,3
rot_g_E(a,b,i)=Qrot_g_E_ini(i,i,fiber_dir(a))
end do
rot_g_E(a,b,4)=Qrot_g_E_ini(1,2,fiber_dir(a))+Qrot_g_E_ini(2,1,fiber_dir(a))
rot_g_E(a,b,5)=Qrot_g_E_ini(2,3,fiber_dir(a))+Qrot_g_E_ini(3,2,fiber_dir(a))
rot_g_E(a,b,6)=Qrot_g_E_ini(3,1,fiber_dir(a))+Qrot_g_E_ini(1,3,fiber_dir(a))


end do


!------------------------------------------------------------------回転終了-------------------------------------

end do


!ミーゼスひずみ(マクロひずみテンソルの1行分に格納)
do a=1,g_econst
do b=1,g_ipn
rot_g_E(a,b,7) = (5.0d-1 * ((rot_g_E(a,b,1) - rot_g_E(a,b,2)) ** 2.0d0 + (rot_g_E(a,b,2) - rot_g_E(a,b,3)) ** 2.0d0 + (rot_g_E(a,b,3) - rot_g_E(a,b,1)) ** 2.0d0 + 6.0d0 * (rot_g_E(a,b,4) ** 2.0d0 + rot_g_E(a,b,5) ** 2.0d0 + rot_g_E(a,b,6) ** 2.0d0))) **5.0d-1
end do
end do




!-------------------------------------------出力------------------------------------------------------
open(10,file='output_rot_g_E.dat')
do a=1,g_econst
write(10,'("element",I5)')a
do b=1,g_ipn
write(10,'("ipn",I5)')b
write(10,'(7F15.7)')(rot_g_E(a,b,j),j=1,7)
end do
end do

close(10)




deallocate(temp2)

write(*,*)"rotation_g_E "


end subroutine