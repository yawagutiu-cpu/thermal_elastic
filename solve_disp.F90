subroutine solve_disp

use parameter1
implicit none
!大出さんのもので確認したけど変更なし
integer i,j,k,l,m,n
  integer column,row,nrhs
allocate(g_disp(g_dim*g_nconst))
allocate(g_node_after(g_dim*g_nconst))

g_disp=0.0d0
g_node_after=0.0d0

do i=1,g_dim*g_nconst
g_disp(i)=deltaT*BTgamma(i)      !ここでのg_dispは変位ではなく力。。。？ 
end do


 open(10,file='output_sl_global_disp0.dat')
    do i=1,g_nconst
      write(10,'(I4,3F25.7)') i, g_disp(i*3-2), g_disp(i*3-1), g_disp(i*3)
    enddo
  close(10)
 

column=3*g_nconst
row=column
nrhs=1
info=0


 call dgetrs(trans, row, nrhs, g_Kmat_LU, column, ipiv_K_g, g_disp, column, info)    !熱弾性による変位
 

!不正な値の確認
print*,info

!変形後の節点座標
do i=1,g_nconst
do j=1,g_dim
g_node_after(3*i-3+j)=g_node(i,j)+g_disp(3*i-3+j)*50   !*g_mag   !ここもあってる！
end do
end do


!---------------------------------------------------------------------------

  !構造物の変位の出力
 open(10,file='output_global_disp.dat')
    do i=1,g_nconst
      write(10,'(I6,3F25.7)') i, g_disp(i*3-2), g_disp(i*3-1), g_disp(i*3)
    enddo
  close(10)
  
  !変形後の構造物の座標の出力
  open(20,file='output_global_node_after.dat')
    do i=1,g_nconst
      write(20,'(I6,3F25.7)') i, g_node_after(3*i-2), g_node_after(3*i-1), g_node_after(3*i)
    enddo
  close(20)


write(*,*) 'FEM Disp COMPLETE!!'

endsubroutine
