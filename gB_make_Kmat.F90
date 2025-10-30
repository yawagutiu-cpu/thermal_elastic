subroutine gB_make_Kmat
!global座標系におけるKmatを作成 ここもコピペでいいはずなんだよなあ、巨視的弾性剛性テンソルの話だから　とりあえずコピペ
use parameter1
implicit none


integer i,j,k,a,b,i1,i2,j1,j2

!配列サイズを指定
allocate(g_tBDmat(g_econst,g_ipn,g_ipn*3,6))
allocate(g_Kmat(g_nconst*3,g_nconst*3))
allocate(g_eKmat(g_econst,g_ipn*3,g_ipn*3))
allocate(g_Kmat_dummy(g_nconst*3,g_nconst*3))

!配列初期化
g_tBDmat=0.0d0
g_Kmat=0.0d0
g_eKmat=0.0d0
g_Kmat=0.0d0
g_Kmat_dummy=0.0d0

!make_eKmat

!tBDの作成
do a=1,g_econst
do b=1,g_ipn
do i=1,g_ipn*3
do j=1,6
do k=1,6
g_tBDmat(a,b,i,j)=g_tBDmat(a,b,i,j)+g_Bmat(a,b,k,i)*g_Amat(a,k,j)*g_det_J(a,b)  !Dmat(擾乱を含めない弾性剛性テンソル)をg_Amat(グローバル座標に合わせた擾乱を含めた弾性剛性テンソル)に変えて計算 detJもここでかける
end do
end do
end do
end do
end do


do a=1,g_econst
do b=1,g_ipn
!eKmatの作成　tBDBもここで作成
do i=1,g_ipn*3
do j=1,g_ipn*3
do k=1,6
g_eKmat(a,i,j)=g_eKmat(a,i,j)+g_tBDmat(a,b,i,k)*g_Bmat(a,b,k,j)
end do
end do
end do
end do
end do


!add_Kmat
!節点ごとに導出した要素剛性マトリックスeKを足し合せて全体剛性マトリックスを作成
do a=1,g_econst
do i1=1,g_ipn !積分点番号の行
do i2=1,g_ipn !積分点番号の列
do j1=1,3 !積分点番号の次元
do j2=1,3 !積分点番号の次元
g_Kmat((g_element(a,i1)-1)*3+j1,(g_element(a,i2)-1)*3+j2)=g_Kmat((g_element(a,i1)-1)*3+j1,(g_element(a,i2)-1)*3+j2)+g_eKmat(a,(i1-1)*3+j1,(i2-1)*3+j2)
end do
end do
end do
end do
end do

g_Kmat_dummy=g_Kmat  !境界条件を適用する前のg_Kmatを保存しておく(solve_globalで使用)

!----------------------データ出力----------------------------------
!if(option_file(11) == 1) then
!    open(221, file = 'Output_g_Kmat_CE.dat', status = 'replace')
!        write(221, '(A)') 'Output_Kmat_CE'
!        write(fmt, '("("I0"F40.25)")') 3*g_nconst
!        do i=1, g_nconst*3
!            write(221, fmt) (g_Kmat(i,j), j=1, g_nconst*3)
!        end do
!    close(221)
!
!    print *, 'made g_Kmatrix (output file : "Output_g_Kmat.dat")'
!
!else
!    print *, 'made g_Kmatrix'
!
!end if

!if(option_file(12) == 1) then
!    open(223, file = 'Output_g_eKmat.dat', status = 'replace')
!        write(223, '(A)') 'Output_eKmat'
!        do a=1, g_econst
!            write(223, '("element", I3)') a
!            do b=1, g_ipn*3
!                write(223, '(24F40.25)') (g_eKmat(a,b,i), i=1, g_ipn*3)
!            end do
!        end do
!    close(223)
!end if

deallocate (g_eKmat)

!------------------------------------------------------------------

call gBA_BoundaryCondition

!make_g_Kmat_LU
allocate(g_Kmat_LU(3*g_nconst,3*g_nconst))
allocate(ipiv_K_g(3*g_nconst))



info = 0

do i=1,3*g_nconst
do j=1,3*g_nconst
g_Kmat_LU(i,j)=g_Kmat(i,j)
end do
end do


call dgetrf(3*g_nconst,3*g_nconst,g_Kmat_LU,3*g_nconst,ipiv_K_g,info)

!不正な値の確認
print*,info
!----------------------データ出力----------------------------------
!if(option_file(21) == 1) then
!    open(21, file = 'Output_g_Kmat_LU.dat', status = 'replace')
!        write(21, '(A)') 'Output_g_Kmat_LU'
!        write(fmt, '("("I0"F20.10)")') g_nconst*3
!        do i=1, g_nconst*3
!            write(21, fmt) (g_Kmat_LU(i,j), j=1, g_nconst*3)
!        end do
!    close(21)
!
!    print *, 'made global Kmatrix_LU (output file : "Output_g_Kmat_LU.dat")'
!
!else
!    print *, 'made global Kmatrix_LU'
!
!end if
!open(11,file='Output_ipiv_K_g.dat',status='replace')
!do i=1,3*g_nconst
!write(11,'(I4,3F20.10)')ipiv_K_g(i)
!end do
!close(11)


deallocate(g_Kmat)
end subroutine


