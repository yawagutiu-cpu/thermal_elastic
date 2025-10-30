subroutine lB_make_Bmat

!ミクロモデルの変位-歪みマトリックスを作成 rotation_local_structureの方針で変更あり
!make_Hmatrix, make_Jmatrix, make_invJmatrix, make_cartesian, make_Bmatrixを統合

use parameter1
implicit none
integer i, j, k, n, m

!配列サイズを指定
allocate (Hmat(ipn,3,ipn))
allocate (Jmat(econst,ipn,dim,dim))
allocate (invJmat(econst,ipn,dim,dim))
allocate (det_J(econst,ipn))
allocate (dNdxyz(econst,ipn,dim,ipn))
allocate (Bmat(econst,ipn,6,ipn*3))

!配列初期化
Hmat    = 0.0d0
Jmat    = 0.0d0
invJmat = 0.0d0
det_J   = 0.0d0
dNdxyz  = 0.0d0
Bmat    = 0.0d0

if(dim==3)then
do m=1,econst
do n=1,ipn

!make_Hmatrix
!形状関数Ni = (1+ξ*ξi)*(1+η*ηi)*(1+ζ*ζi)/8をξηζで偏微分
do i=1,ipn
Hmat(n,1,i)= polar(i,1) * (1.0d0 + gauss(n,2) * polar(i,2)) * (1.0d0 + gauss(n,3) * polar(i,3)) / 8.0d0
Hmat(n,2,i)= (1.0d0 + gauss(n,1) * polar(i,1)) * polar(i,2) * (1.0d0 + gauss(n,3) * polar(i,3)) / 8.0d0
Hmat(n,3,i)= (1.0d0 + gauss(n,1) * polar(i,1)) * (1.0d0 + gauss(n,2) * polar(i,2)) * polar(i,3) / 8.0d0
end do

!make_Jmat
!Jmat作成
do i=1,dim
do j=1,dim
do k=1,ipn
 Jmat(m,n,i,j) = Jmat(m,n,i,j) + Hmat(n,i,k) * node(element(m,k),j)
end do
end do
end do

det_J(m,n) = Jmat(m,n,1,1) * Jmat(m,n,2,2) * Jmat(m,n,3,3) + Jmat(m,n,1,2) * Jmat(m,n,2,3) * Jmat(m,n,3,1) + Jmat(m,n,1,3) * Jmat(m,n,2,1) * Jmat(m,n,3,2) &
                        &- Jmat(m,n,1,3) * Jmat(m,n,2,2) * Jmat(m,n,3,1) - Jmat(m,n,1,1) * Jmat(m,n,2,3) * Jmat(m,n,3,2) - Jmat(m,n,1,2) * Jmat(m,n,2,1) * Jmat(m,n,3,3)
!make_invJmat
          invJmat(m,n,1,1) = (-1.0d0)**(1.0d0 + 1.0d0) * ((Jmat(m,n,2,2) * Jmat(m,n,3,3)) - (Jmat(m,n,2,3) * Jmat(m,n,3,2)))
            invJmat(m,n,1,2) = (-1.0d0)**(1.0d0 + 2.0d0) * ((Jmat(m,n,1,2) * Jmat(m,n,3,3)) - (Jmat(m,n,3,2) * Jmat(m,n,1,3)))
            invJmat(m,n,1,3) = (-1.0d0)**(1.0d0 + 3.0d0) * ((Jmat(m,n,1,2) * Jmat(m,n,2,3)) - (Jmat(m,n,2,2) * Jmat(m,n,1,3)))
            invJmat(m,n,2,1) = (-1.0d0)**(2.0d0 + 1.0d0) * ((Jmat(m,n,2,1) * Jmat(m,n,3,3)) - (Jmat(m,n,3,1) * Jmat(m,n,2,3)))
            invJmat(m,n,2,2) = (-1.0d0)**(2.0d0 + 2.0d0) * ((Jmat(m,n,1,1) * Jmat(m,n,3,3)) - (Jmat(m,n,3,1) * Jmat(m,n,1,3)))
            invJmat(m,n,2,3) = (-1.0d0)**(2.0d0 + 3.0d0) * ((Jmat(m,n,1,1) * Jmat(m,n,2,3)) - (Jmat(m,n,2,1) * Jmat(m,n,1,3)))
            invJmat(m,n,3,1) = (-1.0d0)**(3.0d0 + 1.0d0) * ((Jmat(m,n,2,1) * Jmat(m,n,3,2)) - (Jmat(m,n,3,1) * Jmat(m,n,2,2)))
            invJmat(m,n,3,2) = (-1.0d0)**(3.0d0 + 1.0d0) * ((Jmat(m,n,1,1) * Jmat(m,n,3,2)) - (Jmat(m,n,3,1) * Jmat(m,n,1,2)))
            invJmat(m,n,3,3) = (-1.0d0)**(3.0d0 + 1.0d0) * ((Jmat(m,n,1,1) * Jmat(m,n,2,2)) - (Jmat(m,n,2,1) * Jmat(m,n,1,2)))
do i=1,3
do j=1,3
 invJmat(m,n,i,j) = invJmat(m,n,i,j) / det_J(m,n)
end do
end do

!make_cartesian(dNdxyz)
do i=1,3
do j=1,ipn
do k=1,dim
dNdxyz(m,n,i,j) = dNdxyz(m,n,i,j) + invJmat(m,n,i,k) * Hmat(n,k,j)
end do
end do
end do

!make_Bmatrix
do i=1,ipn
Bmat(m,n,1,3*i-2)=dNdxyz(m,n,1,i)
Bmat(m,n,4,3*i-1) = dNdxyz(m,n,1,i)
Bmat(m,n,6,3*i)  = dNdxyz(m,n,1,i)
Bmat(m,n,4,3*i-2) = dNdxyz(m,n,2,i)
Bmat(m,n,2,3*i-1) = dNdxyz(m,n,2,i)
Bmat(m,n,5,3*i)   = dNdxyz(m,n,2,i)
Bmat(m,n,6,3*i-2) = dNdxyz(m,n,3,i)
Bmat(m,n,5,3*i-1) = dNdxyz(m,n,3,i)
Bmat(m,n,3,3*i)   = dNdxyz(m,n,3,i)
end do
end do
end do

else if(dim==2)then
do m=1,econst
do n=1,ipn
!make Hmatrix
do i=1,ipn
Hmat(n,1,i) = polar(i,1) * (1.0d0 + gauss(n,2) * polar(i,2)) / 4.0d0
Hmat(n,2,i) = polar(i,2) * (1.0d0 + gauss(n,1) * polar(i,1)) / 4.0d0
end do

!make Jmatrix
do i=1,dim
do k=1,ipn
Jmat(m,n,i,1) = Jmat(m,n,i,1) + Hmat(n,i,k) * node(element(m,k),1)
Jmat(m,n,i,2) = Jmat(m,n,i,2) + Hmat(n,i,k) * node(element(m,k),2)   !繊維方向を回転させてないのでこれで正しいはず
end do
end do

det_J(m,n) = Jmat(m,n,1,1) * Jmat(m,n,2,2) - Jmat(m,n,1,2) * Jmat(m,n,2,1)

!make invJmatrix
 invJmat(m,n,1,1) = Jmat(m,n,2,2) / det_J(m,n)
 invJmat(m,n,2,1) = -Jmat(m,n,2,1) / det_J(m,n)
 invJmat(m,n,1,2) = -Jmat(m,n,1,2) / det_J(m,n)
 invJmat(m,n,2,2) = Jmat(m,n,1,1) / det_J(m,n)

do i=1,dim
do j=1,ipn
do k=1,dim
dNdxyz(m,n,i,j) = dNdxyz(m,n,i,j) + invJmat(m,n,i,k) * Hmat(n,k,j) 
end do
end do
end do

!make_Bmatrix
!Two and half dimension   !ここはrotation local structureに関係するので変更の必要あり2023/12/12
do i=1,ipn
Bmat(m,n,1,3*i-2) = dNdxyz(m,n,1,i)
Bmat(m,n,2,3*i-1)   = dNdxyz(m,n,2,i)
Bmat(m,n,4,3*i-1) = dNdxyz(m,n,1,i)
Bmat(m,n,4,3*i-2) = dNdxyz(m,n,2,i)
Bmat(m,n,5,3*i)   = dNdxyz(m,n,2,i)
Bmat(m,n,6,3*i) = dNdxyz(m,n,1,i)
end do
end do
end do
end if

deallocate (Hmat)
deallocate (dNdxyz)

!----------------------データ出力----------------------------------
if(option_file(28) == 1) then
    open(121, file = 'Output_Bmat_CE_2Dhalf.dat', status = 'replace')
        write(121, '(A)') 'output_Bmat_CE_2Dhalf'
        do m=1, econst
            write(121, '("element No. ", I3)') m
            do n=1, ipn
                write(121, '("ipn No. ", I2)') n
                do i=1, 6
                    write(121, '(33F40.25)') (Bmat(m,n,i,j), j=1, 3*ipn)
                end do
            end do
        end do
    close(121)

    print *, 'made Bmatrix_CE (2D half) (output file : "output_Bmat_CE_2Dhalf.dat")'

else
    print *, 'made Bmatrix_CE (2D half)'

end if

end subroutine