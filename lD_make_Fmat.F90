subroutine lD_make_Fmat
!(均質化過程)Fマトリックスを作成　そのままコピペ
use parameter1
implicit none
integer i, i1, j, j1, n, m

!mode == 1, 2 (適合要素, 非適合要素)の場合

!配列サイズ指定
allocate(eFmat(econst,ipn*3,6))
allocate(Fmat(3*nconst,6))

!配列初期化
eFmat=0.0d0
Fmat=0.0d0

!make eFmat [eF]=[B]t*[D]*|J|*w （w=1）
do m=1,econst
do n=1,ipn
do i=1,ipn*3
do j=1,6
eFmat(m,i,j)=eFmat(m,i,j)+tBDmat(m,n,i,j)  !detJはtBDmatに含まれているため不要
end do
end do
end do
end do

!eFmatを足し合わせてFmatを作成
do m=1,econst
do j=1,6
do i1=1,ipn
do j1=1,3
Fmat((element(m,i1)-1)*3+j1, j) = Fmat((element(m,i1)-1)*3+j1, j) + eFmat(m,(i1-1)*3+j1, j)
end do
end do
end do
end do

Fmat=-Fmat 

deallocate(eFmat)

!----------------------データ出力----------------------------------
if(option_file(4) == 1) then
    open(141, file = 'Output_Fmat.dat', status = 'replace')
        write(141, '(A)') 'Output_Fmat'
        do i=1, nconst*3
            write(141, '(6F40.25)') (Fmat(i,j), j=1, 6)
        end do
    close(141)

    print *, 'made Fmatrix (output file : "Output_Fmat.dat")'

else
    print *, 'made Fmatrix'

end if



end subroutine