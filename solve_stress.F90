subroutine solve_stress

use parameter1
implicit none


integer a,b,i,j,k,l,m,n

allocate(l_S(g_econst,g_ipn,econst,ipn,7))
allocate(l_Se(g_econst,g_ipn,econst,7))
allocate(g_S(g_econst,g_ipn,7))
allocate(g_Se(g_econst,7))
allocate(rot_g_E(g_econst,g_ipn,7))


!変数初期化
g_S=0.0d0
l_S=0.0d0
g_Se=0.0d0
l_Se=0.0d0




!マクロ応力・ひずみ

!Σij(y,T) = <aijkl>*Ekl - ΔT*<γij>    巨視的構成式
do a=1,g_econst
do b=1,g_ipn
do k=1,6
do l=1,6
g_S(a,b,k)=g_S(a,b,k)+g_Amat(a,k,l)*g_E(a,b,l)    !<aijkl>*Ekl   周りの熱弾性によって与えられる巨視的ひずみによって生じる応力
end do
end do
end do
end do

!周りの熱弾性だけの応力の出力--------------------------------------------
!要素ごとg_S

allocate(g_Se1(g_econst,6))
g_Se1=0.0d0

do m=1,g_econst
do n=1,g_ipn
do i=1,6
g_Se1(m,i)=g_Se1(m,i)+g_S(m,n,i)/g_ipn
end do
end do
end do

open(10,file='Output_g_Se1.dat')
do i=1,g_econst
    write(10,'(I7,9F15.7)') i, ( g_Se1(i,j), j=1,6)
  end do
close(10)
!-------------------------------------------------------------------


do a=1,g_econst
do b=1,g_ipn
do k=1,6
g_S(a,b,k)=g_S(a,b,k)-deltaT*g_rot_Gamma(a,k)  ! - ΔT*<γij>  ユニット自体の熱弾性によって与えられる応力
end do
end do
end do


!ユニットの熱弾性だけの応力の出力-----------------------------------------
allocate(g_S1(g_econst,g_ipn,6))
allocate(g_Se2(g_econst,6))
g_S1=0.0d0
g_Se2=0.0d0

do  a=1,g_econst
do b=1,g_ipn
do k=1,6
g_S1(a,b,k)=g_S1(a,b,k)-deltaT*g_rot_Gamma(a,k)  ! - ΔT*<γij>  ユニット自体の熱弾性によって与えられる応力
end do
end do
end do

do m=1,g_econst
do n=1,g_ipn
do i=1,6
g_Se2(m,i)=g_Se2(m,i)+g_S1(m,n,i)/g_ipn
end do
end do
end do

open(10,file='Output_g_Se2.dat')
do i=1,g_econst
    write(10,'(I7,9F15.7)') i, ( g_Se2(i,j), j=1,6)
  end do
close(10)
!----------------------------------------------------------------------



!ミーゼス応力
do m=1, g_econst
do n=1, g_ipn
g_S(m,n,7) = (5.0d-1 * ((g_S(m,n,1) - g_S(m,n,2)) ** 2.0d0 + (g_S(m,n,2) - g_S(m,n,3)) ** 2.0d0 + (g_S(m,n,3) - g_S(m,n,1)) ** 2.0d0 + 6.0d0 * (g_S(m,n,4) ** 2.0d0 + g_S(m,n,5) ** 2.0d0 + g_S(m,n,6) ** 2.0d0))) **5.0d-1
end do
end do




!要素ごとg_S
do m=1,g_econst
do n=1,g_ipn
do i=1,6
g_Se(m,i)=g_Se(m,i)+g_S(m,n,i)/g_ipn
end do
end do
end do

!ミーゼス
do m=1,g_econst
g_Se(m,7)=(0.5*((g_Se(m,1)-g_Se(m,2))**2+(g_Se(m,2)-g_Se(m,3))**2+(g_Se(m,3)-g_Se(m,1))**2+6.0*(g_Se(m,4)**2+g_Se(m,5)**2+g_Se(m,6)**2)))**0.5
end do



call rotation_g_Se_th  !積層座標系でのマクロ応力　見たいよなあ？！



call rotation_g_E   !新規！でも、ここでミクロ座標系にかけるマクロひずみの方向を調整する
!ミクロ応力
do a=1,g_econst
do b=1,g_ipn

!微視的応力
do m=1,econst
do n=1,ipn
do i=1,6
do j=1,6
l_S(a,b,m,n,i)=l_S(a,b,m,n,i)+eAmat(m,n,i,j)*rot_g_E(a,b,j)
end do
end do
end do
end do

do m=1,econst
do n=1,ipn
do i=1,6
l_S(a,b,m,n,i)=l_S(a,b,m,n,i)-deltaT*n_gamma(m,n,i)
end do
end do
end do





!ミーゼス応力
do m=1,econst
do n=1,ipn
l_S(a,b,m,n,7)= (0.5 * ((l_S(a,b,m,n,1) - l_S(a,b,m,n,2)) ** 2 + (l_S(a,b,m,n,2) - l_S(a,b,m,n,3)) ** 2 + (l_S(a,b,m,n,3) - l_S(a,b,m,n,1)) ** 2 &
                             + 6.0 * (l_S(a,b,m,n,4) ** 2 + l_S(a,b,m,n,5) ** 2 + l_S(a,b,m,n,6) **2 ))) ** 0.5
end do
end do

!要素ごとの微視的応力
do m=1,econst
do n=1,ipn
do i=1,6
l_Se(a,b,m,i)=l_Se(a,b,m,i)+l_S(a,b,m,n,i)/ipn
end do
end do
end do


!要素ごとのミーゼス応力
do m=1,econst
l_Se(a,b,m,7) = (5.0d-1 * ((l_Se(a,b,m,1) - l_Se(a,b,m,2)) ** 2.0d0 + (l_Se(a,b,m,2) - l_Se(a,b,m,3)) ** 2.0d0 + (l_Se(a,b,m,3) - l_Se(a,b,m,1)) ** 2.0d0 + 6.0d0 * (l_Se(a,b,m,4) ** 2.0d0 + l_Se(a,b,m,5) ** 2.0d0 + l_Se(a,b,m,6) ** 2.0d0))) **5.0d-1
end do

end do
end do

!-----------------データ出力？？？------------------------------------------
filename="Output_g_S.dat"
open(10,file=filename)
  do a=1,g_econst
    write(10,*)'global element number',a
    do b=1,g_ipn
      write(10,*)' global ipn number'
      write(10,'(I7,7F15.7)') b,(g_S(a,b,i), i=1,7)
    end do
  end do
close(10)

open(20,file='output_g_Se.dat')
  do i=1,g_econst
    write(20,'(I7,9F15.7)') i, ( g_Se(i,j), j=1,7)
  end do
close(20)


!ひとまず l_Seを出力しようか  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!出力するセルを選べるようにする←できた2023/12/22
do j=1,cupling_number
write(CWORKe,"(I6)") option_output_Ep(j,1)
write(CWORKi,"(I2)") option_output_Ep(j,2)
filename="output_l_Se(e="//CWORKe(2:6)//",i="//CWORKi(2:2)//").dat"  !ここのCWORKのルールがよくわからん　2023/12/22
open(30,file=filename)
do m=1,econst
write(30,'(I5,7F20.7)')m,(l_Se(option_output_Ep(j,1),option_output_Ep(j,2),m,i),i=1,7)
end do
close(30)
end do

!g_Ee出力
 open(20,file='output_g_Ee.dat')
  do i=1,g_econst
    write(20,'(I7,9F15.7)') i, ( g_Ee(i,j), j=1,7)
  end do
close(20)




!---------------------------------まだ足りない気がする？？？----------------
write(*,*)'FEM Stress COMPLETE!!'
end subroutine

