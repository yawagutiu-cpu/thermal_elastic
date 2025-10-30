subroutine lJ_make_Psi
!îMíeê´ÇÃì¡ê´ä÷êîÉµ
use parameter1
implicit none


integer i,j,k,l,m,n

 allocate(Psi_vec(3*nconst))
  Psi_vec = 0.0d0

do i=1,3*nconst
Psi_vec(i)=Mmat(i)
end do

call dgetrs(trans, 3*nconst, 1, Kmat_LU, 3*nconst, ipiv_K, Psi_vec, 3*nconst, info)

deallocate(Mmat)

 !èoóÕ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(filename,'("Output_Psivector.dat")') 
    open(10,file=filename)
     do m=1,nconst*3
       write(10,'(1F20.7)') Psi_vec(m)
     end do
    close(10)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "Psivec COMPLETE!!"

end subroutine
