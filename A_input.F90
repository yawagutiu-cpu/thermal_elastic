subroutine A_input
!2023/12/16 変更完了？一部弾粘塑性やってた頃のparameterを残しておく　使わんけどな
use parameter1
implicit none

integer i,j,k

!基本情報を読み込み
open(1, file = 'input_basic_info.txt', form = 'formatted')
    read(1,*) name ;    read(1,*) ga
    read(1,*) name ;    read(1,*) pi

    read(1,*) name ;    read(1,*) dim  
    read(1,*) name ;    read(1,*) econst    !モデルに応じて変更
    read(1,*) name ;    read(1,*) nconst    !モデルに応じて変更
             
    read(1,*) name ;    read(1,*) vol  
    read(1,*) name ;    read(1,*) kind

    read(1,*) name ;    read(1,*) g_dim   
    read(1,*) name ;    read(1,*) g_econst        !モデルに応じて変更
    read(1,*) name ;    read(1,*) g_nconst        !モデルに応じて変更
    read(1,*) name ;    read(1,*) g_ply           !新規追加　積層数　モデルに応じて変更
    read(1,*) name ;    read(1,*) g_lambda  
              
    read(1,*) name ;    read(1,*) C_penalty_flag
    read(1,*) name ;    read(1,*) Y_penalty_flag

!積分点の計算
ipn=4*(dim-1)
g_ipn=4*(g_dim-1)

call AA_make_Matrix

e_rate=0.0d0
young=0.0d0
poisson=0.0d0
alpha=0.0d0


    read(1,*) name ;    read(1,*) trigger
    read(1,*) name ;    read(1,*) deg_f(1)   !0deg層　値は90   ミクロのモデルが90deg層が初期方向のため
    read(1,*) name ;    read(1,*) deg_f(2)   !90deg層　値は0
    read(1,*) name ;    read(1,*) deg_f(3)   !45deg層　値は45
    read(1,*) name ;    read(1,*) deg_f(4)   !-45deg層 値は135
    read(1,*) name ;    read(1,*) dt
    read(1,*) name ;    read(1,*) e_rate(1)
    read(1,*) name ;    read(1,*) l_mag
    read(1,*) name ;    read(1,*) g_mag
    read(1,*) name ;    read(1,*) g_end
    read(1,*) name ;    read(1,*) Umat_speed
    
    read(1,*) name ;    read(1,*) deltaT

close(1)

print '("local_dim=",I2)',dim
print '("local_econst=",I4)',econst
print '("local_nconst=",I4)',nconst
print '("global_dim=",I2)',g_dim
print '("global_econst=",I5)',g_econst
print '("global_nconst=",I5)',g_nconst
print '("deg_f(1)=",F6.2)',deg_f(1)
print '("deg_f(2)=",F6.2)',deg_f(2)
print '("deg_f(3)=",F6.2)',deg_f(3)
print '("deg_f(4)=",F6.2)',deg_f(4)
print '("trigger=",I4)',trigger
print '("global_Umat_speed=",F10.4)',Umat_speed
print '("global_end=",F10.4)',g_end
print '("deltaT=",F10.4)',deltaT


e_rate(2)=e_rate(1)

!時間依存均質化の変位関数の方向
dir=dim+dim*(dim-1)/2

!ラジアン変換
do i=1,4
deg_f(i)=deg_f(i)*pi/1.800d2
end do


!local情報の入力
!local節点情報
open(2,file='input_local_node.txt',form='formatted')
 do i=1, nconst
        read(2,*) k, (node(i,j), j=1, 3)
    end do
close(2)

!local要素情報
open(3, file = 'input_local_element.txt', form = 'formatted')
    do i=1, econst
        read(3,*) k, (element(i,j), j=1, ipn)
    end do
close(3)

!localの材料割り当て情報
open(7, file = 'input_local_material_number.txt', form = 'formatted')
    do i=1, econst
        read(7,*) k, material(i)
    end do
close(7)



!global情報の入力
!global節点情報
open(4, file = 'input_global_point.txt', form = 'formatted')
    do i=1, g_nconst
        read(4,*) k, (g_node(i,j), j=1,g_dim)
    end do
close(4)

!global要素情報
open(5, file = 'input_global_element.txt', form = 'formatted')
    do i=1, g_econst
        read(5,*) k, (g_element(i,j), j=1, g_ipn)
    end do
close(5)

!global積層情報 どの要素が何層目に属するかを示す
open(15,file='input_global_pid.txt',form='formatted')
do i=1,g_econst
read(15,*)k,pid(i)
end do
close(15)

!global積層方向情報　　何番めの積層が何度層なのかを示す　PIDと併せてfiber_dirを求める
open(16,file='input_fiber_deg.txt',form='formatted')
do i=1,g_ply       !積層数に応じて変更
read(16,*)k,ply_deg(i)
end do
close(16)

do i=1,g_ply
ply_deg(i)=ply_deg(i)+1.0d0   !積層方向ごとの割り当て番号がfileは0,1,2,3 programは1,2,3,4なので調整
end do

!要素ごとの積層方向計算
do i=1,g_econst
fiber_dir(i)=ply_deg(pid(i))
end do


open(20,file='Output_fiber_dir.dat',status='replace')    !fiber_dirを確認
do i=1,g_econst
write(20,'(I5,1F6.3)')i,(fiber_dir(i))
end do
close(20)
!要素の積層方向は deg_f(fiber_dir(a))　となる　ここでaはglobalの要素番号

open(30,file='Output_deg_f.dat',status='replace')
do i=1,g_econst
write(30,'(I5,1F6.3)')i,(deg_f(fiber_dir(i)))
end do
close(30)




!材料情報をmaterial_propertyから読み込み
   open(6,file='input_material_property.txt',form='formatted')  
    read(6,*) name
    do i=1, kind
        read(6,*) k, young(i), poisson(i), alpha(i,1)
    end do   
    read(6,*) name
    read(6,*) name, ELL
    read(6,*) name, ETT
    read(6,*) name, VLL
    read(6,*) name, VTT
    read(6,*) name, GLT
    read(6,*) name, nonl
    read(6,*) name, hf_p1
    read(6,*) name, hf_p2
    read(6,*) name, hf_p3
!カーボンの線膨張係数
    read(6,*) name, alpha11
    read(6,*) name, alpha33
    
      do i=1,kind  !エポキシは等方性なので熱膨張の係数もすべて同じ
    alpha(i,2) = alpha(i,1)
    alpha(i,3) = alpha(i,1)
  enddo
 !---------------------------------------------------以下、あとで整理して書き換える 　書き換え済み2023/12/16

do i=1,3
do j=1,3
alpha(kind+i,j)=alpha11
end do
end do

do i=1,3
alpha(kind+i,i)=alpha33
end do
 !--------------------------------------------------------------------------
close(6)


!ファイル出力情報を"file_option.txt"から読込み(どのデータを出力するか選べる)
open(9, file = 'input_file_option.txt', form = 'formatted')
    do i=1, 40
        read(9,*) name ;    read(9,*) option_file(i)
    end do
close(9)

!出力する要素・積分点情報を"output_option.txt"から読込み                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!此処の精査する！2023/11/02/16:16
open(10, file = 'input_output_option.txt', form = 'formatted')
    read(10,*) name ;       read(10,*) cupling_number
    allocate (option_output_Ep(cupling_number, 2))
    read(10,*) name
    do i=1, cupling_number
        read(10,*) option_output_Ep(i,1), option_output_Ep(i,2)
    end do
close(10)


!!出力ステップの情報を"step_option.txt"から読込み
!open(11, file = 'input_step_option.txt', form = 'formatted')
!    read(11,*) name ; read(11,*) step_number
!    allocate (option_output_step(step_number))
!    do i=1, step_number
!        read(11,*) k, option_output_step(i)
!    end do
!close(11)

print *, "Completed input step (01)"

end subroutine