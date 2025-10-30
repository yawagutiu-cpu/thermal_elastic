subroutine g0_rotation_Amat
!より簡単に！   角度の正負のみちょっと不安 2023/12/16
use parameter1

implicit none
!3軸方向を向いているユニットセルの均質化したAmatを、
!①積層方向を考慮して回転
!②ダブテール形状を考慮して回転


integer a,b, e, f, i, j, k, l, m, n, p, q,ii,jj
integer i1, i2, i3, i4, j1, j2, j3, j4
integer ip(6), jp(6), kp(6), lp(6)
integer iq(3,3)     !Amatの座標変換用

!配列サイズを指定
allocate(tensor(3,3,3,3))
allocate(rot_tensor(3,3,3,3))
allocate(g_Qmat(3,3))
allocate(rot_Amat(g_econst,g_ipn,6,6))
allocate(g_Qmat_2(3,3))
allocate(g_Amat(g_econst,6,6))
allocate(g_invAmat(g_econst,6,6))





  !Amatの座標変換 というか4階テンソルに割り当てる
   ip(1: 6) = (/ 1, 2, 3, 1, 2, 3 /)
   jp(1: 6) = (/ 1, 2, 3, 2, 3, 1 /)
   kp(1: 6) = (/ 1, 2, 3, 1, 2, 3 /)
   lp(1: 6) = (/ 1, 2, 3, 2, 3, 1 /)
   
   iq(1, 1: 3) = (/1, 4, 6/)
   iq(2, 1: 3) = (/4, 2, 5/)
   iq(3, 1: 3) = (/6, 5, 3/)
   






do a=1,g_econst
do b=1,g_ipn
tensor=0.0d0
do i=1,3
do j=1,3
do k=1,3
do l=1,3
tensor(i,j,k,l)=Amat_temp(iq(i,j),iq(k,l))
                !tensor(1,1,1,1) = Amat_temp(1,1)
                !tensor(2,2,1,1) = Amat_temp(2,1)
                !tensor(3,3,1,1) = Amat_temp(3,1)
                !tensor(1,1,2,2) = Amat_temp(1,2)
                !tensor(2,2,2,2) = Amat_temp(2,2)
                !tensor(3,3,2,2) = Amat_temp(3,2)
                !tensor(1,1,3,3) = Amat_temp(1,3)
                !tensor(2,2,3,3) = Amat_temp(2,3)
                !tensor(3,3,3,3) = Amat_temp(3,3)

                !tensor(2,1,2,1) = Amat_temp(4,4)
                !tensor(1,2,2,1) = Amat_temp(4,4)
                !tensor(2,1,1,2) = Amat_temp(4,4)
                !tensor(1,2,1,2) = Amat_temp(4,4)

                !tensor(3,2,3,2) = Amat_temp(5,5)
                !tensor(2,3,3,2) = Amat_temp(5,5)
                !tensor(3,2,2,3) = Amat_temp(5,5)
                !tensor(2,3,2,3) = Amat_temp(5,5)

                !tensor(3,1,3,1) = Amat_temp(6,6)
                !tensor(1,3,3,1) = Amat_temp(6,6)
                !tensor(3,1,1,3) = Amat_temp(6,6)
                !tensor(1,3,3,1) = Amat_temp(6,6)
end do
end do
end do
end do

!①積層方向を考慮して1軸周りの回転（積層方向はdeg_f(fiber_dir(a)))
g_Qmat=0.0d0

g_Qmat(1,1)=1
g_Qmat(2,2)=cos(deg_f(fiber_dir(a)))
g_Qmat(2,3)=-sin(deg_f(fiber_dir(a)))
g_Qmat(3,2)=sin(deg_f(fiber_dir(a)))
g_Qmat(3,3)=cos(deg_f(fiber_dir(a)))

do j1=1,3
do j2=1,3
do j3=1,3
do j4=1,3
rot_tensor(j1,j2,j3,j4)=0.0d0
do i1=1,3
do i2=1,3
do i3=1,3
do i4=1,3
rot_tensor(j1,j2,j3,j4)=rot_tensor(j1,j2,j3,j4)+tensor(i1,i2,i3,i4)*g_Qmat(i1,j1)*g_Qmat(i2,j2)*g_Qmat(i3,j3)*g_Qmat(i4,j4)
end do
end do
end do
end do
end do
end do
end do
end do


  !ここまでのテンソル回転分をtensorに上書き
     do j1=1, 3
       do j2=1, 3
         do j3=1, 3
           do j4=1, 3
             tensor(j1, j2, j3, j4) = rot_tensor(j1, j2, j3, j4)
           end do
         end do
       end do
     end do
     
     
!②要素の傾きに応じた3軸周りの回転
 g_Qmat_2=0.0d0
 g_Qmat_2(1,1)=cos(-th(a))
 g_Qmat_2(1,2)=-sin(-th(a))
 g_Qmat_2(2,1)=sin(-th(a))
 g_Qmat_2(2,2)=cos(-th(a))
 g_Qmat_2(3,3)=1
 
 do j1=1,3
 do j2=1,3
 do j3=1,3
 do j4=1,3
rot_tensor(j1,j2,j3,j4)=0.0d0
do i1=1,3
do i2=1,3
do i3=1,3
do i4=1,3
rot_tensor(j1,j2,j3,j4)=rot_tensor(j1,j2,j3,j4)+tensor(i1,i2,i3,i4)*g_Qmat_2(i1,j1)*g_Qmat_2(i2,j2)*g_Qmat_2(i3,j3)*g_Qmat_2(i4,j4)
end do
end do
end do
end do
end do
end do
end do
end do


!2階のAmatに格納
  do jj=1, 6
       do ii=1, 6
         i = ip(ii)
         j = jp(ii)
         k = kp(jj)
         l = lp(jj)
         rot_Amat(a,b,ii, jj) = (rot_tensor(i, j, k, l) + rot_tensor(i, j, l, k)) / 2.0 
       end do
     end do


     
     !global用のAmatに格納
     do i=1, 6
       do j=1, 6
         g_Amat(a, i, j) = rot_Amat(a,b,i, j)      !global座標にそろったAmatであるg_Amatができた
          
       end do
     end do




end do
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!回転終了!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 open(10, file='output_g_Amat.dat')
     do a=1, g_econst
       write(10, '("element", I5)') a
       do i=1, 6
         write(10, '(6F20.10)') (g_Amat(a, i, j), j=1, 6)
       end do
     end do
   close(10)
   
   open(20, file='output_theata.dat')
     do a=1, g_econst
       write(20, '(I5, 1F20.10)') a, th(a)
     end do
   close(20)
   print *, 'g_Amat'


  !Amat_inv---------------------------
   allocate(Imat(6, 6))
   allocate(Amat_dummy(6, 6))
   allocate(Amat_inv_dummy(6, 6))
   
   Imat = 0.0d0
   
   do i=1, 6
     Imat(i, i) = 1.0d0
   end do
   
   do a=1, g_econst
     Amat_dummy = 0.0d0
     Amat_inv_dummy = 0.0d0
     do i=1, 6
       do j=1, 6 
         Amat_dummy(i, j) = g_Amat(a, i, j)
       end do
     end do
     
     
     call make_inverse(Amat_dummy, Amat_inv_dummy, 6)
     
     do i=1, 6
       do j=1, 6
         g_invAmat(a, i, j) = Amat_inv_dummy(i, j)
       end do
     end do
   end do
   
   
   deallocate(tensor)
   deallocate(rot_tensor)
   deallocate(rot_Amat)
   deallocate(Imat)
   deallocate(Amat_dummy)
   deallocate(Amat_inv_dummy)
   
   
   print *, '/   file      g_Amat, theate'
   print *, 'rotation_Amat complete'
   
!   read *
   
   end subroutine


