subroutine C_rotation_local_structure
!現状使っていない　2023/12/16

use parameter1
implicit none

integer i,j,k,n

!配列サイズを指定
allocate(node_dummy(nconst,3))
allocate(alpha_dummy(3))            !新規追加
allocate(node_f(nconst,3))
allocate(rot_mat(4,3,3))
allocate(node_rot(3,nconst,3))
allocate(rot_mat_rev(3,3,3))

!配列初期化
node_dummy = 0.0d0
alpha_dummy=0.0d0   !新規追加
rot_mat = 0.0d0
node_rot = 0.0d0


!ローカル節点nodeをx1方向に座標変換
do n=1,nconst
node_dummy(n,1)=node(n,3)
node_dummy(n,2)=-node(n,2)
node_dummy(n,3)=node(n,1)
end do

node=0.0d0

do n=1,nconst
do i=1,3
node(n,i)=node_dummy(n,i)
end do
end do


!alphaも変換

!異方性材料のalphaをx1方向に座標変換
alpha_dummy(1) = alpha(kind+1,3)
alpha_dummy(2) = alpha(kind+1,1)
alpha_dummy(3) = alpha(kind+1,2)

do n=1,3
alpha(kind+1,n)=0.0d0
end do

do n=1,3
alpha(kind+1,n)=alpha_dummy(n)
end do

!!deg   これややこしいしいらんくね
!!x1軸
!rot_mat(1,1,1)=cos(deg)
!rot_mat(1,1,2)=sin(deg)
!rot_mat(1,1,3)=0.0d0
!!x2軸
!rot_mat(1,2,1)=sin(deg)
!rot_mat(1,2,2)=cos(deg)
!rot_mat(1,2,3)=0.0d0
!!x3軸
!rot_mat(1,3,1) = 0.0d0
!rot_mat(1,3,2) = 0.0d0
!rot_mat(1,3,3) = 1.0d0
!
!do n=1,nconst
!do i=1,3
!do j=1,3
!node(n,i)=node(n,i)+rot_mat(1,i,j)*node_dummy(n,j)    
!end do
!end do
!end do
!
!!alphaも!
!do i=1,3
!do j=1,3
!alpha(kind+1,i)=alpha(kind+1,i)+rot_mat(1,i,j)*alpha_dummy(j)        !とても不安　あってる？
!end do
!end do


open(31, file = 'output_rotation_local_node.dat', status = 'replace')
    do n=1, nconst
        write(31, '(I5, 3F15.7)') n, (node(n,i), i=1, 3)
    end do
close(31)

!以下、一要素では使わない？
!node_f=0.0d0
!
!do n=1,nconst
!do i=1,3
!node_f(n,i)=node_dummy(n,i)
!end do
!end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!こっから先わかんね！！！！！！！！！！！！後回し？？？？
!
!
!!deg_fごとに回転行列を作成(i=1:0deg, i=2:90deg, i=3:45deg)
!!x1-x2平面上での(x3軸を中心とした)回転行列→nodeをx1方向について座標変換しているので実質x1-x3平面での(x2軸を中心とした)回転行列
!do i=1,3
!rot_mat(i+1,1,1)=cos(deg_f(i))
!rot_mat(i+1,1,2)=-sin(deg_f(i))
!rot_mat(i+1,1,3)=0.0d0
!rot_mat(i+1,2,1)=sin(deg_f(i))
!rot_mat(i+1,2,2)=cos(deg_f(i))
!    rot_mat(i+1,2,3) = 0.0d0
!    rot_mat(i+1,3,1) = 0.0d0
!    rot_mat(i+1,3,2) = 0.0d0
!    rot_mat(i+1,3,3) = 1.0d0
!
!    if(deg_flag(i) == 1) then  !cos(90deg)=1,sin(90deg)=1に補正
!        rot_mat(i,1,1) = 1.0d0
!        rot_mat(i,1,2) = 0.0d0
!        rot_mat(i,2,1) = 0.0d0
!        rot_mat(i,2,2) = 1.0d0
!    end if
!end do
!
!
!
!!deg_fごとに節点を回転させて、node_rotに代入
!do i=1,3
!do n=1,nconst
!do j=1,3
!do k=1,3
! node_rot(i,n,j) = node_rot(i,n,j) + rot_mat(i,k,j) * node_f(n,k)
!end do
!end do
!end do
!end do
!
!!---------deg_fに回転させた節点を逆回転して戻すときに使う用のrot_mat_rev--------------------------------------------------
!!配列初期化
!rot_mat_rev = 0.0d0
!
!do i=1,3
!!rot_mat_revは回転しない場合
!
!    rot_mat_rev(i,1,1) = cos(-deg_f(i))
!    rot_mat_rev(i,1,2) = -sin(-deg_f(i))
!    rot_mat_rev(i,1,3) = 0.0d0
!    rot_mat_rev(i,2,1) = sin(-deg_f(i))
!    rot_mat_rev(i,2,2) = cos(-deg_f(i))
!    rot_mat_rev(i,2,3) = 0.0d0
!    rot_mat_rev(i,3,1) = 0.0d0
!    rot_mat_rev(i,3,2) = 0.0d0
!    rot_mat_rev(i,3,3) = 1.0d0
!    
!    if(deg_flag(i) == 1) then  !cos(90deg)=1,sin(90deg)=1に補正
!        rot_mat_rev(i,1,1) = 1.0d0
!        rot_mat_rev(i,1,2) = 0.0d0
!        rot_mat_rev(i,2,i) = 0.0d0
!        rot_mat_rev(i,2,2) = 1.0d0
!    end if
!end do

deallocate (node_dummy)
deallocate (rot_mat)

end subroutine


!
!subroutine C_rotation_local_structure
!!node,alphaともに1軸方向に回転させる
!use parameter1
!implicit none
!
!integer i,j,k,l,m,n,a,b
!
!allocate(node_dummy(nconst,3))
!allocate(alpha_dummy(3))            
!allocate(node_f(nconst,3))     !辻川さんのnode_temp
!allocate(rot_mat(4,3,3))     !辻川さんのQmat_3rot_1~4dir
!allocate(node_rot(4,nconst,3))   !辻川さんのnode1~4
!allocate(rot_mat_rev(3,3,3))
!!allocate(node_temp(nconst,3))
!
!allocate(alpha_rot(4,4,3))     !新規追加　回転させた後のalpha 一応出しとく　alphaって　　nodeと違ってもしかして使わないんじゃないかと思ったり思わなかったり　マクロで回転させたGammaは必要そうだけど
!
!
!!配列初期化
!node_dummy = 0.0d0
!alpha_dummy=0.0d0   !新規追加
!rot_mat = 0.0d0
!node_rot = 0.0d0
!
!!ローカル節点nodeをx1方向に座標変換
!do n=1,nconst
!node_dummy(n,1)=node(n,3)
!node_dummy(n,2)=-node(n,2)
!node_dummy(n,3)=node(n,1)
!end do
!
!node=0.0d0
!
!do n=1,nconst
!do i=1,3
!node(n,i)=node_dummy(n,i)
!end do
!end do
!
!
!!alphaも変換
!
!!異方性材料のalphaをx1方向に座標変換
!alpha_dummy(1) = alpha(kind+1,3)
!alpha_dummy(2) = alpha(kind+1,1)
!alpha_dummy(3) = alpha(kind+1,2)
!
!do n=1,3
!alpha(kind+1,n)=0.0d0
!end do
!
!do n=1,3
!alpha(kind+1,n)=alpha_dummy(n)
!end do
!
!
!
!!!90deg積層　　3軸方向が繊維方向なのでそのままでOK--------------------------------------------
!!do i=1,nconst
!!do j=1,3
!!node_rot(2,i,j)=node(i,j)
!!end do
!!end do
!!
!!!alphaも回転
!!do i=1,3
!!alpha_rot(2,4,i)=alpha(4,i)
!!end do
!!
!!!0deg積層   3軸方向が繊維方向であるユニットセルを2軸方向を繊維方向に向き変える---------------
!!do i=1,nconst
!!node_rot(1,i,1)=node(i,1)
!!node_rot(1,i,2)=node(i,3)
!!node_rot(1,i,3)=node(i,2)
!!end do
!!
!!alpha_rot(1,4,1)=alpha(4,1)
!!alpha_rot(1,4,2)=alpha(4,3)
!!alpha_rot(1,4,3)=alpha(4,2)
!!
!!
!!
!!!+-45deg積層------------------------------------------------多分これであってる気がする　まあ検証用0,90だからはまだ使わないけど
!!!1軸中心に45deg回転させる
!!do n=1, nconst  !判りやすいようにとりあえず0degに向ける
!!    node_dummy(n,1) = node(n,1)
!!    node_dummy(n,2) = node(n,3)
!!    node_dummy(n,3) = node(n,2)
!!end do
!!
!!node_f = 0.0d0
!!
!!do n=1, nconst
!!    do i=1, 3
!!        node_f(n,i) = node_dummy(n,i)
!!    end do
!!end do
!!
!!
!!do i=3,4    !2軸繊維方向のユニットセルをglobalの1軸中心にdeg_fだけ回転させる　　　ここが不安
!!    rot_mat(i,1,1) =1.0d0
!!    rot_mat(i,1,2) =0.0d0
!!    rot_mat(i,1,3) = 0.0d0
!!    rot_mat(i,2,1) =0.0d0
!!    rot_mat(i,2,2) =cos(deg_f(i))
!!    rot_mat(i,2,3) =-sin(deg_f(i))
!!    rot_mat(i,3,1) = 0.0d0
!!    rot_mat(i,3,2) =sin(deg_f(i))
!!    rot_mat(i,3,3) =cos(deg_f(i))
!!end do
!!
!!
!!!deg_fごとに節点を回転させて，node_rotに代入
!!do i=3,4
!!    do n=1, nconst
!!        do j=1, 3
!!            do k=1, 3
!!                node_rot(i,n,j) = node_rot(i,n,j) + rot_mat(i,k,j) * node_f(n,k)
!!            end do
!!        end do
!!    end do
!!end do
!!
!!!alpha
!!alpha_dummy(1)=alpha(4,1)
!!alpha_dummy(2)=alpha(4,3)
!!alpha_dummy(3)=alpha(4,2)
!!
!!do i=3,4
!!do j=1,3
!!do k=1,3
!!alpha_rot(i,4,j)=alpha_rot(i,4,j)+rot_mat(i,k,j)*alpha_dummy(k)
!!end do
!!end do
!!end do
!!
!!
!!
!!
!!!!---------deg_fに回転させた節点を逆回転して戻すときに使う用のrot_mat_rev　　ここまだ精査できてない--------------------------------------------------
!!!!配列初期化
!!!rot_mat_rev = 0.0d0
!!!
!!!do i=1, 3
!!!    !rot_mat_rev(1,i,j)は回転しない場合
!!!    rot_mat_rev(i,1,1) = cos(-deg_f(i))
!!!    rot_mat_rev(i,1,2) = -sin(-deg_f(i))
!!!    rot_mat_rev(i,1,3) = 0.0d0
!!!    rot_mat_rev(i,2,1) = sin(-deg_F(i))
!!!    rot_mat_rev(i,2,2) = cos(-deg_f(i))
!!!    rot_mat_rev(i,2,3) = 0.0d0
!!!    rot_mat_rev(i,3,1) = 0.0d0
!!!    rot_mat_rev(i,3,2) = 0.0d0
!!!    rot_mat_rev(i,3,3) = 1.0d0
!!!    
!!!    if(deg_flag(i) == 1) then  !cos(90deg)=1,sin(90deg)=1に補正
!!!        rot_mat_rev(i,1,1) = 1.0d0
!!!        rot_mat_rev(i,1,2) = 0.0d0
!!!        rot_mat_rev(i,2,i) = 0.0d0
!!!        rot_mat_rev(i,2,2) = 1.0d0
!!!    end if
!!!end do
!!
!!
!!!!今後の諸々の計算のためにnodeを1方向に向ける
!!! do n=1,nconst
!!!        node_temp(n,1) = node(n,3)
!!!        node_temp(n,2) = node(n,1)
!!!        node_temp(n,3) = node(n,2)
!!!    end do
!!!
!!!  node = 0.0d0                !nodeを座標変換に対応するようにnode_tempに保存したので初期化
!!!    
!!!    do n=1,nconst
!!!    do i=1,dim
!!!    node(n,i)=node_temp(n,i)
!!!    end do
!!!    end do
!!    
!!
!!deallocate (node_dummy)
!!deallocate (rot_mat)
!
!end subroutine
!
