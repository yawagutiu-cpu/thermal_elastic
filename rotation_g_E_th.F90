subroutine rotation_g_E_th
!要素の傾き座標でのひずみを見たい！thだけ逆回転させる
use parameter1
implicit none


integer a,b,i,j,k,l,p,q
integer i1,i2,j1,j2,k1,k2
real(8) temp(3,3)                             !計算結果を一時的に保存しておくため
allocate(g_E_mat(3,3))
allocate(rot_g_E_ini(3,3))

allocate(rot_g_E_th(g_econst,g_ipn,7))
allocate(rot_g_Ee_th(g_econst,7))

rot_g_E_th=0.0d0


do a=1,g_econst

 

!要素の傾きに応じた3軸周りの回転  逆方向回転なので-thでなくth
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





!格納
do i=1,3
rot_g_E_th(a,b,i)=rot_g_E_ini(i,i)
end do
rot_g_E_th(a,b,4)=rot_g_E_ini(1,2)+rot_g_E_ini(2,1)
rot_g_E_th(a,b,5)=rot_g_E_ini(2,3)+rot_g_E_ini(3,2)
rot_g_E_th(a,b,6)=rot_g_E_ini(3,1)+rot_g_E_ini(1,3)


end do



end do


!ミーゼスひずみ(マクロひずみテンソルの1行分に格納)
do a=1,g_econst
do b=1,g_ipn
rot_g_E_th(a,b,7) = (5.0d-1 * ((rot_g_E_th(a,b,1) - rot_g_E_th(a,b,2)) ** 2.0d0 + (rot_g_E_th(a,b,2) - rot_g_E_th(a,b,3)) ** 2.0d0 + (rot_g_E_th(a,b,3) - rot_g_E_th(a,b,1)) ** 2.0d0 + 6.0d0 * (rot_g_E_th(a,b,4) ** 2.0d0 + rot_g_E_th(a,b,5) ** 2.0d0 + rot_g_E_th(a,b,6) ** 2.0d0))) **5.0d-1
end do
end do

rot_g_Ee_th=0.0d0

!要素ごとのマクロひずみ
do a=1,g_econst
do b=1,g_ipn
do i=1,6
rot_g_Ee_th(a,i)=rot_g_Ee_th(a,i)+rot_g_E_th(a,b,i)/g_ipn
end do
end do
end do

!要素ごとのミーゼスひずみ
do a=1, g_econst
rot_g_Ee_th(a,7) = (5.0d-1 * ((rot_g_Ee_th(a,1) - rot_g_Ee_th(a,2)) ** 2.0d0 + (rot_g_Ee_th(a,2) - rot_g_Ee_th(a,3)) ** 2.0d0 + (rot_g_Ee_th(a,3) - rot_g_Ee_th(a,1)) ** 2.0d0 + 6.0d0 * (rot_g_Ee_th(a,4) ** 2.0d0 + rot_g_Ee_th(a,5) ** 2.0d0 + rot_g_Ee_th(a,6) ** 2.0d0))) **5.0d-1
end do
!-------------------------------------------出力------------------------------------------------------

!rotation_g_Ee_th出力
  open(20,file='Output_rot_g_Ee_th.dat',status='replace')
        do a=1, g_econst
            write(20,'(I7,7F15.7)')a ,(rot_g_Ee_th(a,i), i=1,7)
        end do
    close(20)



deallocate(g_E_mat)
deallocate(rot_g_E_ini)
deallocate(rot_g_E_th)
deallocate(rot_g_Ee_th)

write(*,*)"rotation_g_E_th "


end subroutine
