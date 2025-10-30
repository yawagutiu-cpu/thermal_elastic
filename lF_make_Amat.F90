subroutine lF_make_Amat
!均質化弾性合成テンソルを作成 そのままコピペ
use parameter1
implicit none
integer i,j,k,m,n

allocate(BXmat(econst,ipn,6,6))
allocate(eAmat(econst,ipn,6,6))
allocate(Amat_temp(6,6))
allocate(invAmat_temp(6,6))
allocate(Imat(6,6))

!配列初期化
eAmat=0.0d0
BXmat=0.0d0
Amat_temp=0.0d0
invAmat_temp=0.0d0
Imat=0.0d0

!単位行列の作成
do i=1,6
Imat(i,i)=1.0d0
end do

!BXmatrixの作成 BX=[B][eX]
do m=1,econst
do n=1,ipn
do i=1,6
do j=1,6
do k=1,ipn*3
BXmat(m,n,i,j)=BXmat(m,n,i,j)+Bmat(m,n,i,k)*eXmat(m,k,j)
end do
end do
end do
end do
end do

!eAmatの計算　eAmat=Dmat*(Imat+BXmat)   (2.15)
do m=1,econst
do n=1,ipn
do i=1,6
do j=1,6
do k=1,6
eAmat(m,n,i,j)=eAmat(m,n,i,j)+Dmat(material(m),i,k)*(Imat(k,j)+BXmat(m,n,k,j))
end do
end do
end do
end do
end do

!Amatの計算
do m=1,econst
do n=1,ipn
do i=1,6
do j=1,6
Amat_temp(i,j)=Amat_temp(i,j)+eAmat(m,n,i,j)*det_J(m,n)
end do
end do
end do
end do

Amat_temp=Amat_temp/Vol !体積平均をとる

!本来0であるべき成分の誤差を修正
do i=1, 6
    do j=1, 6
        if(Amat_temp(i,j) < 1.0d-5) Amat_temp(i,j) = 0.0d0
    end do
end do

!ランダム配列の場合に影響してくるせん断成分を0にする
!do i=1, 3
!   do j=1, 3
!       Amat_temp(i, 3+j) = 0.0d0
!       Amat_temp(3+j, i) = 0.0d0
!       Amat_temp(i+j-1, 3+j) = 0.0d0
!       Amat_temp(3+i, i+j-1) = 0.0d0
!   end do
!end do

!invAmatの算出
call make_inverse(Amat_temp,invAmat_temp,6)

deallocate (eXmat)
deallocate (BXmat)
deallocate (Imat)  

!----------------------データ出力----------------------------------
if(option_file(8) == 1) then
    open(161, file = 'Output_Amat.dat', status = 'replace')
        write(161, '(A)') 'Output_Amat'
        do i=1, 6
            write(161, '(6F50.35)') (Amat_temp(i,j), j=1, 6)
        end do
    close(161)

    print *, 'made Amatrix (output file : "Output_Amat.dat")'

else
    print *, 'made Amatrix'

end if

if(option_file(9) == 1) then
    open(162, file = 'Output_invAmat.dat', status = 'replace')
        write(162, '(A)') 'output_invAmat'
        do i=1, 6
            write(162, '(6F40.25)') (invAmat_temp(i,j), j=1, 6)
        end do
    close(162)
end if

open(163, file = 'Output_eAmat.dat', status = 'replace')
    write(163, '(A)') 'Output_eAmat'
    do m=1, econst
        write(163, '("element"I5)') m
        do n=1, ipn
            write(163, '("ipn"I3)') n
            do i=1, 6
                write(163, '(6F50.35)') (eAmat(m,n,i,j), j=1, 6)
            end do
        end do
    end do
close(163)

end subroutine