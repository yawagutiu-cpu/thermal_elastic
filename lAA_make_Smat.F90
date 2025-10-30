subroutine lAA_make_Smat

!炭素繊維3軸方向それぞれについてコンプライアンス行列を作成 
use parameter1
implicit none

integer i,j,k

!配列サイズを指定

allocate (Smat(6,6))
allocate (invSmat(3,6,6))
allocate (invSmat_p(6,6))


!配列初期化
invSmat=0.0d0

!Smat(コンプライアンス行列)の作成 今回は繊維方向が3軸なので「炭素繊維3軸」のみ用いる　
!炭素繊維1軸
Smat=0.0d0
Smat(1,1)=1.0d0/ELL
Smat(1,2)=-VLL/ELL
Smat(1,3)=Smat(1,2)
Smat(2,1)=Smat(1,2)
Smat(2,2)=1.0d0/ETT
Smat(2,3)=-VTT/ETT
Smat(3,1)=Smat(1,2)
Smat(3,2)=Smat(2,3)
Smat(3,3)=Smat(2,2)

Smat(4,4)=1.0d0/GLT
Smat(5,5)=2.0d0*(Smat(2,2)-Smat(2,3))
Smat(6,6)=Smat(4,4)

invSmat_p=0.0d0
call make_inverse(Smat,invSmat_p,6)

do i=1,6
do j=1,6
invSmat(1,i,j)=invSmat_p(i,j) 
end do
end do

!炭素繊維2軸
Smat = 0.0d0
Smat(1,1) = 1.0d0 / ETT
Smat(1,2) = -VLL / ELL
Smat(1,3) = -VTT / ETT
Smat(2,1) = Smat(1,2)
Smat(2,2) = 1.0d0 / ELL
Smat(2,3) = Smat(1,2)
Smat(3,1) = Smat(1,3)
Smat(3,2) = Smat(1,2)
Smat(3,3) = Smat(1,1)

Smat(4,4) = 1.0d0 / GLT
Smat(5,5) = Smat(4,4)
Smat(6,6) = 2.0d0 * (Smat(3,3) - Smat(3,1))

invSmat_p = 0.0d0
call make_inverse(Smat,invSmat_p,6)

do i=1, 6
    do j=1, 6
        invSmat(2,i,j) = invSmat_p(i,j)
    end do
end do

!炭素繊維が3軸方向の時の物性
Smat = 0.0d0
Smat(1,1) = 1.0d0 / ETT
Smat(1,2) = -VTT / ETT
Smat(1,3) = -VLL / ELL
Smat(2,1) = Smat(1,2)
Smat(2,2) = Smat(1,1)
Smat(2,3) = Smat(1,3)
Smat(3,1) = Smat(1,3)
Smat(3,2) = Smat(2,3)
Smat(3,3) = 1.0d0 / ELL

Smat(4,4) = 2.0d0 * (Smat(1,1) - Smat(1,2))
Smat(5,5) = 1.0d0 / GLT
Smat(6,6) = Smat(5,5)

invSmat_p = 0.0d0
call make_inverse(Smat,invSmat_p,6)

do i=1, 6
    do j=1, 6
        invSmat(3,i,j) = invSmat_p(i,j)
    end do
end do

deallocate (Smat)
deallocate (invSmat_p)
        
print *, 'made Smatrix (output file : "invSmat.dat")'
        
                 
end subroutine