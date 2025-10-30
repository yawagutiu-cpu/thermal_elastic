subroutine make_BTgamma
!変位を求めるための布石？？   大出さんのものを参考に一要素の時とは書き換えた！一要素の時は基本泉崎さんを参考にしていたけど間違ってるかもって先生にいわれたから
use parameter1

implicit none

integer a,b,i,j,k,l,m,n
!allocate(e_BTgamma(g_econst,g_ipn,3*g_ipn))
allocate(e_BTgamma(g_econst,3*g_ipn))
allocate(BTgamma(3*g_nconst))

!macro_stress=gamma*deltaT
e_BTgamma=0.0d0
BTgamma=0.0d0

do a=1,g_econst
do b=1,g_ipn
do i=1,3*g_ipn
do j=1,6
e_BTgamma(a,i)=e_BTgamma(a,i)+g_Bmat(a,b,j,i)*g_rot_Gamma(a,j)*g_det_J(a,b)   !要素ごとのBTgamma detJかけるのなんで？？← [B]の転置×[Gamma] makeBTR参照
!赤座さんのBTRを参考に、e_BTGammaを2回のテンソルにしたらそれっぽい値が出た。BTGammaの値は大きいのにg_dispは小さくなった。なんで？2023/12/13
end do
end do
end do
end do

do a=1,g_econst
do b=1,g_ipn
do i=1,3
BTgamma(3*g_element(a,b)+i-3)=BTgamma(3*g_element(a,b)+i-3)+e_BTgamma(a,3*b+i-3)  !ここはあってそう! 2023/12/13
end do
end do
end do

!---------------------------------------出力---------------------------------------------------------
open(10,file='output_BTGamma.dat')
  do i=1,g_nconst
    write(10,'(I5,3F15.7)') i, BTgamma(i*3-2), BTgamma(i*3-1), BTgamma(i*3)
  end do
close(10)

open(20,file='output_e_BTgamma.dat')
  do a=1,g_econst
    write(20,*) 'element number', a
    do b=1,g_ipn
!      write(20,*) '  ipn number', b
!      do i=1,g_ipn
        write(20,'(I4,3F15.7)') b, e_BTgamma(a,3*b-2), e_BTgamma(a,3*b-1), e_BTgamma(a,3*b)
      enddo
    enddo
!  enddo
close(20)
deallocate(e_BTgamma)

write(*,*) "BTGamma COMPLETE!!"

end subroutine

