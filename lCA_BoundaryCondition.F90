subroutine lCA_BoundaryCondition

!ローカルモデルの境界条件(C対称性、Y周期性)を適用 そのままコピペ

use parameter1
implicit none

integer i,j,k

if(C_penalty_flag==1)then
 !C対称性(繊維の軸中心に関して点対象)
   open(1311, file = 'input_local_C-penalty.txt', form = 'formatted')
        read(1311, *) name
        read(1311, *) C1
        read(1311, *) C2
        read(1311, *) C3
        read(1311, *) C4

        read(1311, *) name
        read(1311, *) pair

allocate(pn_C(pair,2))

       read(1311,*) name
        read(1311,*) lambda
pn_C=0.0d0

do i=1,pair
          read(1311,*) k, pn_C(i,1), pn_C(i,2)
end do

close(1311)

  print *, 'C1 = ',C1
    print *, 'C2 = ',C2
    print *, 'C3 = ',C3
    print *, 'C4 = ',C4
    print *, 'lambda = ', lambda


!C1, C2, C3, C4を固定
do i=1,3

       Kmat(3*C1+i-3, 3*C1+i-3) = Kmat(3*C1+i-3, 3*C1+i-3) + lambda
        Kmat(3*C2+i-3, 3*C2+i-3) = Kmat(3*C2+i-3, 3*C2+i-3) + lambda
        Kmat(3*C3+i-3, 3*C3+i-3) = Kmat(3*C3+i-3, 3*C3+i-3) + lambda
        Kmat(3*C4+i-3, 3*C4+i-3) = Kmat(3*C4+i-3, 3*C4+i-3) + lambda
end do

!ペアリング
do i=1,pair
do j=1,3
            Kmat(3*(pn_C(i,1)-1)+j, 3*(pn_C(i,1)-1)+j) = Kmat(3*(pn_C(i,1)-1)+j, 3*(pn_c(i,1)-1)+j) + lambda
            Kmat(3*(pn_C(i,2)-1)+j, 3*(pn_C(i,1)-1)+j) = Kmat(3*(pn_C(i,2)-1)+j, 3*(pn_C(i,1)-1)+j) + lambda
            Kmat(3*(pn_C(i,1)-1)+j, 3*(pn_C(i,2)-1)+j) = Kmat(3*(pn_C(i,1)-1)+j, 3*(pn_C(i,2)-1)+j) + lambda
            Kmat(3*(pn_C(i,2)-1)+j, 3*(pn_C(i,2)-1)+j) = Kmat(3*(pn_C(i,2)-1)+j, 3*(pn_C(i,2)-1)+j) + lambda    
        end do
    end do
    
   deallocate(pn_C)

    print *, "C_penalty BC is given"

end if

if(Y_penalty_flag==1)then   !ローカルモデルがハーフユニットででない場合
    !Y周期性(ユニットセルごとに周期構造が並ぶ)
open(1312,file = 'input_Y-penalty.txt', form = 'formatted')
      read(1312, *) name
        read(1312, *) pair

allocate(pn_Y(pair,2))

   read(1312,*) name
        read(1312,*) lambda

do i=1,pair
 read(1312, *) k, pn_Y(i, 1), pn_Y(i, 2)
end do
close(1312)

!ペアリング
do i=1,pair
do j=1,3
Kmat(3*(pn_Y(i,1)-1)+j,3*(pn_Y(i,1)-1)+j)=Kmat(3*(pn_Y(i,1)-1)+j,3*(pn_Y(i,1)-1)+j)+lambda
Kmat(3*(pn_Y(i,2)-1)+j,3*(pn_Y(i,1)-1)+j)=Kmat(3*(pn_Y(i,2)-1)+j,3*(pn_Y(i,1)-1)+j)-lambda
Kmat(3*(pn_Y(i,1)-1)+j,3*(pn_Y(i,2)-1)+j)=Kmat(3*(pn_Y(i,1)-1)+j,3*(pn_Y(i,2)-1)+j)-lambda
Kmat(3*(pn_Y(i,2)-1)+j,3*(pn_Y(i,2)-1)+j)=Kmat(3*(pn_Y(i,2)-1)+j,3*(pn_Y(i,2)-1)+j)+lambda
end do
end do

deallocate(pn_Y)

 print *, "Y_penalty is given"


end if

!ペナルティ法を適用したKmatの出力
if(option_file(29) == 1) then
    open(1311,file='Output_Kmat(p).dat')
        write(fmt,'("("I0"10F20.7)")') nconst*3
        do i=1,nconst*3
            write(1311,fmt) (Kmat(i,j), j=1,nconst*3)
        end do
    close(1311)
end if


end subroutine


