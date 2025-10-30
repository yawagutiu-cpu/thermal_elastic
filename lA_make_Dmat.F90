subroutine lA_make_Dmat

!ミクロモデルの応力-ひずみマトリックスの作成 汎用性を高める　改正済　2023/12/16

use parameter1
implicit none

integer i,j,k



!配列サイズ指定
allocate(Dmat(kind+3,6,6))
!初期化
Dmat=0.0d0

!Dmat作成（等方性材料の種類ごと）
do k=1,kind

!左上
do i=1,3
do j=1,3
if(i==j)then
 Dmat(k, i, j) = young(k) * (1.0d0 - poisson(k)) / (1.0d0 + poisson(k)) / (1.0d0 - 2.0d0 * poisson(k))
else
Dmat(k, i, j) = young(k) * poisson(k) / (1.0d0 + poisson(k)) / (1.0d0 - 2.0d0 * poisson(k))
end if
end do
end do

!右下
do i=4,6
Dmat(k, i, i) = young(k) / (1.0d0 + poisson(k)) / 2.0d0
end do

end do

call lAA_make_Smat !異方性材料
!Dmatのうち炭素繊維を適用する要素の成分をコンプライアンス行列(invSmat)で置き換え
do k=1,3
do i=1,6
do j=1,6
Dmat(kind+k, i, j) = invSmat(k, i, j)! 3軸方向なので...Dmat(5)に対応させなきゃダメ2023/12/16
end do
end do
end do

deallocate (invSmat)

!----------------------データ出力----------------------------------
if(option_file(1) == 1) then
    open(111, file = 'output_Dmat.dat', status = 'replace')
        write(111, '(A)') 'output_Dmat'
        do i=1, kind+3
            write(111, '("kind", I3)') i
            do j=1, 6
                write(111, '(6F40.25)') (Dmat(i,j,k), k=1, 6)
            end do
        end do
    close(111)

    print *, 'made Dmatrix (output file : "Dmat.dat")'

else
    print *, 'made Dmatrix'

end if

end subroutine