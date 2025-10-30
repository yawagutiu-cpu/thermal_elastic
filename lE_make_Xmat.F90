subroutine lE_make_Xmat
!特性関数χを作成
use parameter1
implicit none
integer i,j,m,n

!配列サイズを指定
allocate(eXmat(econst,ipn*3,6))
allocate(Xmat(nconst*3,6))
allocate(Xmat_temp(nconst*3,6))

!配列初期化
Xmat=0.0d0
eXmat=0.0d0
Xmat_temp=0.0d0

Xmat_temp=Fmat

!Xmatの計算　[k][x]=[f]
call dgetrs(trans,nconst*3,6,Kmat_LU,nconst*3,ipiv_K,Xmat_temp,nconst*3,info)
!call solve_equation(Kmat, Fmat, Xmat, nconst*3, nconst*3, 6)   !solve_equationを使うとoverflowした

Xmat=Xmat_temp

!eXmatの計算
do m=1,econst
do n=1,ipn
do i=1,3
do j=1,6
eXmat(m, 3*(n-1)+i, j) = Xmat(3*(element(m, n)-1)+i, j)
end do
end do
end do
end do

deallocate(Fmat)
deallocate(Xmat_temp)

!----------------------データ出力----------------------------------
if(option_file(5) == 1) then
    open(151, file = 'Output_Xmat.dat', status = 'replace')
        write(151, '(A)') 'Output_Xmat'
        do m=1, econst
            write(151, '("element", I3)') m
            do i=1, ipn*3
                write(151, '(6F15.7)') (eXmat(m,i,j), j=1, 6)
            end do
        end do
    close(151)

    print *, 'made Xmatrix (output file : "Output_Xmat.dat")'

else
    print *, 'made Xmatrix'

end if


end subroutine

