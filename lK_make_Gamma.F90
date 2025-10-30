subroutine lK_make_Gamma

use parameter1

implicit none

integer i,j,k,l,m,n

 allocate(e_Psi(econst,3*ipn))
  allocate(DB(6,3*ipn))
  allocate(Dalpha(econst,ipn,6))
  allocate(DBPsi(econst,ipn,6))
  allocate(n_gamma(econst,ipn,6))
  allocate(gamma(6))
  allocate(macro_alpha(6))

 e_Psi = 0.0d0
  Dalpha = 0.0d0
  DBPsi = 0.0d0
  n_gamma = 0.0d0
  gamma = 0.0d0
  macro_alpha = 0.0d0

!e_Psi
do m=1,econst
do n=1,ipn
do i=1,3
e_Psi(m,3*n-3+i)=Psi_vec(3*element(m,n)-3+i)
end do
end do
end do

!Gamma_vec
do m=1,econst
do n=1,ipn
DB=0.0d0


do i=1,6
do j=1,3*ipn
do k=1,6
DB(i,j)=DB(i,j)+Dmat(material(m),i,k)*Bmat(m,n,k,j)  !応力_ひずみ
end do
end do
end do

do i=1,6
do j=1,6
Dalpha(m,n,i)= Dalpha(m,n,i) + Dmat(material(m),i,j) * alpha(material(m),j)
end do
do j=1,3*ipn
DBPsi(m,n,i)=DBPsi(m,n,i)+DB(i,j)*e_Psi(m,j)
end do
end do
end do
end do

!積分点Gamma
do m=1,econst
do n=1,ipn
do i=1,6
n_gamma(m,n,i)=Dalpha(m,n,i)-DBPsi(m,n,i)  !熱弾性ツースケールスライド P2右下Gammaの式
end do
end do
end do

!全体gamma
do m=1,econst
do n=1,ipn
do i=1,6
gamma(i)=gamma(i)+n_gamma(m,n,i)*det_J(m,n)
end do
end do
end do


!均質化
do i=1,6
gamma(i)=gamma(i)/vol
end do

do i=1,6
do j=1,6
macro_alpha(i)=macro_alpha(i)+invAmat_temp(i,j)*gamma(j)
end do
end do

  if(option_file(13) == 1) then
    open(10,file='Output_gamma.dat')
     write(FMT,'("("I0"F25.15)")') 6
      do i=1,6
        write(10,FMT) gamma(i)
      end do
    close(10)
   end if
  

   open(10,file='Output_n_gamma.dat')
    write(FMT,'("("I0"F20.10)")') 6
     do m=1,econst
      write(10,*)'eleme=',m
      do n=1,ipn
       write(10,*)'ipn=',n
       do i=1,6
        write(10,FMT) n_gamma(m,n,i)
       end do
      end do
     end do
   close(10)


   open(10,file='Output_macro_alpha.dat')
    do i=1,6 
      write(10,'(I4,F20.10)') i,macro_alpha(i)
    end do
   close(10)

    
  write(*,*)'end gamma'
 
  deallocate(e_Psi)
  deallocate(DB)
  deallocate(DBPsi)
  deallocate(Dalpha)
  
end subroutine
