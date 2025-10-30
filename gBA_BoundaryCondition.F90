subroutine gBA_BoundaryCondition
!penaltyを与える節点とその方向を指定
use parameter1
implicit none

integer i,j,k,l,m,n,a,b
integer i1,i2,j1,j2
 
!配列初期化
g_pd = 0.0d0
x_temp = 0.0d0
y_temp = 0.0d0
z_temp = 0.0d0
x_middle = 0.0d0
y_middle = 0.0d0


!x座標の最大値をx_tempに保存
do i=1, g_nconst
    if(g_node(i,1) > x_temp) then
        x_temp = g_node(i,1)
    end if
end do

!y座標の最大値をy_tempに保存
do i=1, g_nconst
    if(g_node(i,2) > y_temp) then
        y_temp = g_node(i,2)
    end if
end do

!z座標の最大値をz_tempに保存
do i=1, g_nconst
    if(g_node(i,3) > z_temp) then
        z_temp = g_node(i,3)
    end if
end do

x_middle = x_temp / 2.0d0
y_middle = y_temp / 2.0d0

print *, 'x_temp = ', x_temp
print *, 'y_temp = ', y_temp
print *, 'z_temp = ', z_temp
print *, 'x_middle = ', x_middle
print *, 'y_middle = ', y_middle





 !-----------------------------------------------------境界条件入力------------
    !x_spcのinput                                   境界での節点番号   
    !x方向固定する節点を指定　今回はダブテールが半分のモデルなのでx軸に垂直な断面をx方向に固定
    open(90,file='x_spc.txt',form='formatted')
    read(90,*) name
    read(90,*) spc_x_number
    read(90,*) name
    print*,spc_x_number
     allocate(x_spc(spc_x_number))
    do i=1,spc_x_number
        read(90,*) x_spc(i)
    end do
    close(90)
    
    open(90, file='output_x_spc.dat')
    do i=1,spc_x_number
        write(90,'(F20.10)') (x_spc(i))
    end do
    close(90)
!    
        !y_spcのinput                                   境界での節点番号   
    !y方向固定する節点を指定　今回はダブテールが半分のモデルなのでx軸に垂直な断面をx方向に固定
    open(90,file='y_spc.txt',form='formatted')
    read(90,*) name
    read(90,*) spc_y_number
    read(90,*) name
    print*,spc_y_number
     allocate(y_spc(spc_y_number))
    do i=1,spc_y_number
        read(90,*) y_spc(i)
    end do
    close(90)
    
         !z_spcのinput                                   境界での節点番号   
    !z方向固定する節点を指定　今回はダブテールが半分のモデルなのでz軸に垂直な断面をz方向に固定
    open(90,file='z_spc.txt',form='formatted')
    read(90,*) name
    read(90,*) spc_z_number
    read(90,*) name
    print*,spc_z_number
     allocate(z_spc(spc_z_number))
    do i=1,spc_z_number
        read(90,*) z_spc(i)
    end do
    close(90)

!全方向固定　とりあえず回転しないように端っこの２点を固定？
 open(90,file='xyz_spc.txt',form='formatted')
    read(90,*) name
    read(90,*) spc_xyz_number
   close(90)
   allocate(xyz_spc(spc_xyz_number))
     open(90,file='xyz_spc.txt',form='formatted')
    read(90,*) name
    read(90,*) spc_xyz_number
    read(90,*) name
    do i=1,spc_xyz_number
        read(90,*) xyz_spc(i)
    end do
    close(90)
  
!--------------------------境界条件付加---------------------------------------------

    !x軸(1軸)に垂直な面を固定
    !x方向固定   x固定なので"*3-2"する
    do i=1,spc_x_number
        g_Kmat(x_spc(i)*3-2,x_spc(i)*3-2) = g_Kmat(x_spc(i)*3-2,x_spc(i)*3-2) + g_lambda
   end do
   print*,'add_x_penalty'
!   
!     !y軸(2軸)に垂直な面を固定
    !y方向固定   y固定なので"*3-1"する
    do i=1,spc_y_number
        g_Kmat(y_spc(i)*3-1,y_spc(i)*3-1) = g_Kmat(y_spc(i)*3-1,y_spc(i)*3-1) + g_lambda
   end do
   print*,'add_y_penalty'
   
     !z軸(1軸)に垂直な面を固定
    !z方向固定   z固定なので"*3"する
    do i=1,spc_z_number
        g_Kmat(z_spc(i)*3,z_spc(i)*3) = g_Kmat(z_spc(i)*3,z_spc(i)*3) + g_lambda
   end do
   print*,'add_z_penalty'
!   
! do i=1,spc_xyz_number    !全方向固定
!        g_Kmat(xyz_spc(i)*3,xyz_spc(i)*3) = g_Kmat(xyz_spc(i)*3,xyz_spc(i)*3) +g_lambda
!        g_Kmat(xyz_spc(i)*3-1,xyz_spc(i)*3-1) = g_Kmat(xyz_spc(i)*3-1,xyz_spc(i)*3-1) + g_lambda
!        g_Kmat(xyz_spc(i)*3-2,xyz_spc(i)*3-2) = g_Kmat(xyz_spc(i)*3-2,xyz_spc(i)*3-2) + g_lambda
!    end do
!   print*,'add_xyz_penalty'



!-------------------------------------------------------------------------------------------------------------------

print*,'add_global_penalty'


end subroutine