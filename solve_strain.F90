subroutine solve_strain


use parameter1
implicit none

integer i,j,k,l,m,n,a,b

allocate(disp_dummy1(g_nconst,g_dim))
allocate(disp_dummy2(g_econst,g_ipn*g_dim))
allocate(g_E(g_econst,g_ipn,7))
allocate(g_Ee(g_econst,7))

disp_dummy1=0.0d0
disp_dummy2=0.0d0
g_E=0.0d0


!1次元配列→2次元配列に入れ替え
do a=1,g_nconst
do b=1,3
disp_dummy1(a,b)=g_disp(a*3-3+b)
end do
end do

!要素ごとに並び変える
do a=1,g_econst
do b=1,g_ipn
do k=1,3
disp_dummy2(a,b*3-3+k)=disp_dummy1(g_element(a,b),k)
end do
end do
end do

!global積分点でのlocalマクロひずみ
do a=1,g_econst
do b=1,g_ipn
do i=1,6
do j=1,3*g_ipn
g_E(a,b,i)=g_E(a,b,i)+g_Bmat(a,b,i,j)*disp_dummy2(a,j)   !Bmatをかけなおすことでひずみに変換 
end do
end do
end do
end do


!ミーゼスひずみ(マクロひずみテンソルの1行分に格納)
do a=1,g_econst
do b=1,g_ipn
g_E(a,b,7) = (5.0d-1 * ((g_E(a,b,1) - g_E(a,b,2)) ** 2.0d0 + (g_E(a,b,2) - g_E(a,b,3)) ** 2.0d0 + (g_E(a,b,3) - g_E(a,b,1)) ** 2.0d0 + 6.0d0 * (g_E(a,b,4) ** 2.0d0 + g_E(a,b,5) ** 2.0d0 + g_E(a,b,6) ** 2.0d0))) **5.0d-1
end do
end do

!g_E出力確認
open(30,file='Output_g_E.dat')
do a=1,g_econst
write(30,*)'g_econst',a
do b=1,g_ipn
write(30,*)'g_ipn'
write(30,'(I7,9F15.7)')b ,(g_E(a,b,i),i=1,7)
end do
end do
close(30)

!----------------------------------------------------此処までは大出さんとおなじ-----------------------------------------------
g_Ee=0.0d0   !初期化！


!要素ごとのマクロひずみ
do a=1,g_econst
do b=1,g_ipn
do i=1,6
g_Ee(a,i)=g_Ee(a,i)+g_E(a,b,i)/g_ipn
end do
end do
end do

!要素ごとのミーゼスひずみ
do a=1, g_econst
g_Ee(a,7) = (5.0d-1 * ((g_Ee(a,1) - g_Ee(a,2)) ** 2.0d0 + (g_Ee(a,2) - g_Ee(a,3)) ** 2.0d0 + (g_Ee(a,3) - g_Ee(a,1)) ** 2.0d0 + 6.0d0 * (g_Ee(a,4) ** 2.0d0 + g_Ee(a,5) ** 2.0d0 + g_Ee(a,6) ** 2.0d0))) **5.0d-1
end do



!g_Ee出力
  open(20,file='Output_g_Ee.dat',status='replace')
        do a=1, g_econst
            write(20,'(I7,7F15.7)')a ,(g_Ee(a,i), i=1,7)
        end do
    close(20)
!------------------------ここから先確認用
call rotation_g_E_th

print*,'FEM strain completed'
end subroutine