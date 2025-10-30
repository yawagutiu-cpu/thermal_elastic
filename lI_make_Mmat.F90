subroutine lI_make_Mmat

use parameter1

implicit none

integer i,j,k,l,m,n

allocate(Mmat(3*nconst))
allocate(e_Mmat(econst,3*ipn))

Mmat=0.0d0
e_Mmat=0.0d0




!eMmat
do m=1,econst
do n=1,ipn
do i=1,3*ipn
do j=1,6
e_Mmat(m,i)=e_Mmat(m,i)+tBDmat(m,n,i,j)*alpha(material(m),j)   !Ç∑Ç≈Ç…detJÇÕÇ©ÇØÇƒÇ†ÇÈÇÃÇ≈äÑà§
end do
end do
end do
end do

!Mmat
do m=1,econst
do n=1,ipn
do k=1,3
Mmat(3*element(m,n)+k-3)=Mmat(3*element(m,n)+k-3)+e_Mmat(m,3*n+k-3)
end do
end do
end do

 !MmatÇÃèoóÕ
   open(10,file='output_Mmat.dat')
    write(FMT,'("("I0"F20.10)")') 6
     do i=1,nconst
       write(10,FMT) Mmat(3*i-2),Mmat(3*i-1),Mmat(3*i)
     end do
   close(10)

 write(*,*)'end Mmat'
  
  deallocate(e_Mmat)
  
  
end subroutine
