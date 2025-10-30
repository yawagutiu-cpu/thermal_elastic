subroutine gA_make_Bmat

!コピペ　特段の違いはない
!グローバル座標におけるBmatを計算
use parameter1
implicit none
integer a,b,i,j,k,m,n

!配列サイズを指定
allocate(g_Hmat(g_ipn,g_dim,g_ipn))
allocate(g_Jmat(g_econst,g_ipn,g_dim,g_dim))
allocate(g_invJmat(g_econst,g_ipn,g_dim,g_dim))
allocate(g_det_J(g_econst,g_ipn))
allocate(g_dNdxyz(g_econst,g_ipn,g_dim,g_ipn))
allocate(g_Bmat(g_econst,g_ipn,6,g_ipn*g_dim))

!配列初期化
g_Hmat = 0.0d0
g_Jmat = 0.0d0
g_invJmat = 0.0d0
g_det_J = 0.0d0
g_dNdxyz = 0.0d0
g_Bmat = 0.0d0

do a=1,g_econst
do b=1,g_ipn
!make_g_Hmat
do i=1,g_ipn
g_Hmat(b,1,i)=polar(i,1)*(1.0d0+gauss(b,2)*polar(i,2))*(1.0d0+gauss(b,3)*polar(i,3))/8.0d0
g_Hmat(b,2,i)=(1.0d0+gauss(b,1)*polar(i,1))*polar(i,2)*(1.0d0+gauss(b,3)*polar(i,3))/8.0d0
g_Hmat(b,3,i)=(1.0d0+gauss(b,1)*polar(i,1))*(1.0d0+gauss(b,2)*polar(i,2))*polar(i,3)/8.0d0
end do

!make_g_Jmat
do i=1,g_dim
do j=1,g_dim
do k=1,g_ipn
g_Jmat(a,b,i,j)=g_Jmat(a,b,i,j)+g_Hmat(b,i,k)*g_node(g_element(a,k),j)
end do
end do
end do





!det_J
      g_det_J(a,b) = g_Jmat(a,b,1,1) * g_Jmat(a,b,2,2) * g_Jmat(a,b,3,3) + g_Jmat(a,b,1,2) * g_Jmat(a,b,2,3) * g_Jmat(a,b,3,1) + g_Jmat(a,b,1,3) * g_Jmat(a,b,2,1) * g_Jmat(a,b,3,2) &
        &- g_Jmat(a,b,1,3) * g_Jmat(a,b,2,2) * g_Jmat(a,b,3,1) - g_Jmat(a,b,1,1) * g_Jmat(a,b,2,3) * g_Jmat(a,b,3,2) - g_Jmat(a,b,1,2) * g_Jmat(a,b,2,1) * g_Jmat(a,b,3,3)

!make_g_invJmat
     g_invJmat(a,b,1,1) = (-1.0d0)**(1.0d0 + 1.0d0) * ((g_Jmat(a,b,2,2) * g_Jmat(a,b,3,3)) - (g_Jmat(a,b,2,3) * g_Jmat(a,b,3,2)))
        g_invJmat(a,b,1,2) = (-1.0d0)**(1.0d0 + 2.0d0) * ((g_Jmat(a,b,1,2) * g_Jmat(a,b,3,3)) - (g_Jmat(a,b,3,2) * g_Jmat(a,b,1,3)))
        g_invJmat(a,b,1,3) = (-1.0d0)**(1.0d0 + 3.0d0) * ((g_Jmat(a,b,1,2) * g_Jmat(a,b,2,3)) - (g_Jmat(a,b,2,2) * g_Jmat(a,b,1,3)))
        g_invJmat(a,b,2,1) = (-1.0d0)**(2.0d0 + 1.0d0) * ((g_Jmat(a,b,2,1) * g_Jmat(a,b,3,3)) - (g_Jmat(a,b,3,1) * g_Jmat(a,b,2,3)))
        g_invJmat(a,b,2,2) = (-1.0d0)**(2.0d0 + 2.0d0) * ((g_Jmat(a,b,1,1) * g_Jmat(a,b,3,3)) - (g_Jmat(a,b,3,1) * g_Jmat(a,b,1,3)))
        g_invJmat(a,b,2,3) = (-1.0d0)**(2.0d0 + 3.0d0) * ((g_Jmat(a,b,1,1) * g_Jmat(a,b,2,3)) - (g_Jmat(a,b,2,1) * g_Jmat(a,b,1,3)))
        g_invJmat(a,b,3,1) = (-1.0d0)**(3.0d0 + 1.0d0) * ((g_Jmat(a,b,2,1) * g_Jmat(a,b,3,2)) - (g_Jmat(a,b,3,1) * g_Jmat(a,b,2,2)))
        g_invJmat(a,b,3,2) = (-1.0d0)**(3.0d0 + 1.0d0) * ((g_Jmat(a,b,3,1) * g_Jmat(a,b,1,2)) - (g_Jmat(a,b,3,2) * g_Jmat(a,b,1,1)))
        g_invJmat(a,b,3,3) = (-1.0d0)**(3.0d0 + 1.0d0) * ((g_Jmat(a,b,1,1) * g_Jmat(a,b,2,2)) - (g_Jmat(a,b,2,1) * g_Jmat(a,b,1,2)))
do i=1,3
do j=1,3
g_invJmat(a,b,i,j)=g_invJmat(a,b,i,j)/g_det_J(a,b)
end do
end do


!make_g_cartesian(dNdxyz)
do i=1,g_ipn
do j=1,g_dim
do k=1,g_dim
g_dNdxyz(a,b,k,i)=g_dNdxyz(a,b,k,i)+g_invJmat(a,b,k,j)*g_Hmat(b,j,i)
end do
end do
end do

!make_g_Bmat
do i=1,g_ipn
g_Bmat(a,b,1,3*i-2)=g_dNdxyz(a,b,1,i)
g_Bmat(a,b,4,3*i-2)=g_dNdxyz(a,b,2,i)
g_Bmat(a,b,6,3*i-2)=g_dNdxyz(a,b,3,i)

g_Bmat(a,b,4,3*i-1)=g_dNdxyz(a,b,1,i)
g_Bmat(a,b,2,3*i-1)=g_dNdxyz(a,b,2,i)
g_Bmat(a,b,5,3*i-1)=g_dNdxyz(a,b,3,i)

g_Bmat(a,b,6,3*i)=g_dNdxyz(a,b,1,i)
g_Bmat(a,b,5,3*i)=g_dNdxyz(a,b,2,i)
g_Bmat(a,b,3,3*i)=g_dNdxyz(a,b,3,i)

end do


end do
end do

deallocate (gauss)
deallocate (polar)
deallocate (g_Hmat)
deallocate (g_Jmat)
deallocate (g_invJmat)
deallocate (g_dNdxyz)

!----------------------データ出力----------------------------------
!if(option_file(17) == 1) then
!    open(211, file = 'output_g_Bmat.dat', status = 'replace')
!        write(211, '(A)') 'output_g_Bmat'
!        do a=1, g_econst
!            write(211, '("global_element", I3)') a
!            do b=1, g_ipn
!                write(211, '("global_ipn", I2)') b
!                do j=1, 6
!                    write(211, '(33F40.25)') (g_Bmat(a,b,j,k), k=1, 3*g_ipn)
!                end do
!            end do
!        end do
!    close(211)
!
!   open(11, file = 'output_g_detJ.dat', status = 'replace')
!        write(11, '(A)') 'output_g_detJ'
!        do a=1, g_econst
!                    write(11, '(I5,33F40.25)')a, (g_det_J(a,k), k=1, g_ipn)
!                end do
!
!    close(11)
!    print *, 'made global Bmatrix (output file : "output_g_Bmat.dat")'
!
!else
!    print *, 'made global Bmatrix'
!
!end if


end subroutine
