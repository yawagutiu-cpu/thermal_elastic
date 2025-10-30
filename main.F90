program main
!ダブテールの熱弾性　辻川さんのsrcを逐次計算に戻して、俺の熱弾性ツースケール一要素と組み合わせる
use parameter1
implicit none
integer m,n,i
call A_input
call B_coordinates
!call C_rotation_local_structure

!ミクロ座標系でのマトリクス計算
call lA_make_Dmat
call lB_make_Bmat
call lC_make_Kmat
call lD_make_Fmat
call lE_make_Xmat
call lF_make_Amat
!熱弾性
call lI_make_Mmat
call lJ_make_Psi
call lK_make_Gamma

!座標変換　此処が難しい！
call make_theata   !要素の傾き角を計算する
call g0_rotation_Amat
call g02_rotation_Gamma

!ユニットセル自体の熱弾性を確認--------------------
allocate(l_S1(econst,ipn,6))
allocate(l_S1e(econst,6))
l_S1=0.0d0
l_S1e=0.0d0

do m=1,econst
do n=1,ipn
do i=1,6
l_S1(m,n,i)=l_S1(m,n,i)-deltaT*n_gamma(m,n,i)
end do
end do
end do

do m=1,econst
do n=1,ipn
do i=1,6
l_S1e(m,i)=l_S1e(m,i)+l_S1(m,n,i)/ipn
end do
end do
end do


open(10,file='l_S1.dat')
do m=1,econst
write(10,'(I5,6F20.7)')m,(l_S1e(m,i),i=1,6)
end do
close(10)
!---------------------------------------------------

!マクロ
call gA_make_Bmat
call gB_make_Kmat

!計算
call solve_thermal_elastic

print*,'Program is completed'


end program