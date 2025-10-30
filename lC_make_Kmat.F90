subroutine lC_make_Kmat
!ローカルモデルの全体剛性マトリックス そのままコピペ

use parameter1
implicit none


integer  i, j, k, n, m, i1, i2, j1, j2

!配列サイズを指定
allocate (tBDmat(econst,ipn,ipn*3,6))
allocate (tBDBmat(econst,ipn*3,ipn*3))
allocate (eKmat(econst,ipn*3,ipn*3))
allocate (Kmat(nconst*3,nconst*3))
allocate (ipiv_K(nconst*3))
allocate (Kmat_LU(nconst*3, nconst*3))

!配列初期化
tBDmat  = 0.0d0
tBDBmat = 0.0d0
eKmat   = 0.0d0
Kmat    = 0.0d0
Kmat_LU = 0.0d0
ipiv_K  = 0.0d0

!make ekmatrix
do m=1,econst
do n=1,ipn
![B]t[D]マトリックスの作成
do i=1,ipn*3
do j=1,6
do k=1,6
tBDmat(m,n,i,j)=tBDmat(m,n,i,j)+Bmat(m,n,k,i)*Dmat(material(m),k,j)*det_J(m,n)
end do
end do
end do
end do
end do

!tBDBmatの作成
do m=1,econst
do n=1,ipn
do i=1,ipn*3
do j=1,ipn*3
do k=1,6
tBDBmat(m,i,j)=tBDBmat(m,i,j)+tBDmat(m,n,i,k)*Bmat(m,n,k,j)
end do
end do
end do
end do
end do

!要素剛性マトリックスの作成
![K] = +[B]t[D][B] * det[J] *　w
do m=1,econst
do i=1,ipn*3
do j=1,ipn*3
eKmat(m,i,j)=eKmat(m,i,j)+tBDBmat(m,i,j)
end do
end do
end do

!add Kmat
do m=1,econst
do i1=1,ipn
do i2=1,ipn
do j1=1,3
do j2=1,3
Kmat((element(m,i1)-1)*3+j1, (element(m,i2)-1)*3+j2 ) = Kmat( (element(m,i1)-1)*3+j1, (element(m,i2)-1)*3+j2 ) + eKmat(m,(i1-1)*3+j1,(i2-1)*3+j2)
end do
end do
end do
end do
end do

!----------------------データ出力----------------------------------
if(option_file(2)==1)then
    open(131, file = 'Output_Kmat_CE.dat', status = 'replace') !適合要素？
        write(131, '(A)') 'Output_Kmat_CE'
        write(fmt,'("("I0"6F15.7)")') nconst*3
        do i=1, nconst*3
            write(131, fmt) (Kmat(i,j), j=1, nconst*3)
        end do
    close(131)

    print *, 'made Kmatrix_CE (output file : "Output_Kmat_CE.dat")'

else
    print *, 'made Kmatrix_CE'

end if

if(option_file(3) == 1) then
    open(133, file = 'Output_eKmat.dat', status = 'replace')
        write(133, '(A)') 'Output_eKmat'
        do m=1, econst
            write(133, '("element", I3)') m
            do n=1, ipn*3
                write(133, '(6F15.7)') (eKmat(m,n,i), i=1, ipn*3)
            end do
        end do
    close(133)
end if



!-----------------------------------------------------------------------------

call lCA_BoundaryCondition !ユニットセル単位での境界条件を指定

!make_Kmat_LU
do i=1,nconst*3
do j=1,nconst*3
Kmat_LU(i,j)=Kmat(i,j)
end do
end do

call dgetrf(nconst*3, nconst*3, Kmat_LU, nconst*3, ipiv_K, info)

!----------------------データ出力----------------------------------
if(option_file(7) == 1) then
    open(132, file = 'Output_Kmat_LU.dat', status = 'replace')
        write(fmt, '("("I0"F70.35)")') nconst*3
        do i=1, nconst*3
            write(132, fmt) (Kmat_LU(i,j), j=1, nconst*3)
        end do
    close(132)

    print *, 'made Kmatrix_LU (output file : "Output_Kmat_LU.dat")'

else
    print *, 'made Kmatrix_LU'

end if

deallocate (Jmat)
deallocate (invJmat)
deallocate (tBDBmat)
deallocate (eKmat)


end subroutine