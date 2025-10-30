subroutine make_theata
!辻川さんのcalculate_th2を参考に作成　　森先輩に確認したい
use parameter1
implicit none

 integer i, j, k, l, a, b, temp
  !  double precision temp2
    
    double precision,allocatable :: node_e1(:)
    double precision,allocatable :: node_e2(:)
  
    allocate (th_ge(g_ipn))   !"parameter1"にinteger,allocatable:: th_ge(:) 追加
    allocate (temp_th(g_econst))    !"parameter1"にdouble precision,allocatable :: temp_th(:) 追加
    allocate (check_th(g_econst,2))
    allocate (node_e1(2))
    allocate (node_e2(2))
    allocate (th_single(g_econst,2))
    allocate(th(g_econst))
  check_th=0.0d0
  
  do a=1,g_econst
  
  !-----------------------------------------------------------------各要素の積分点番号を繊維方向を計算しやすくするために並び替え----------------------------------------------------------------------------------------
    do i=1, g_ipn
      th_ge(i) = i    !th_geに1～8を順番に代入
    end do
    

   !1.z座標の大きいものから順に　th_ge(1),th_ge(2),...と並べ替え
    do k=1, g_ipn-1  !昇順が完全に完了させるためにglobal_ipn-1回繰り返す．1回の操作で隣り合う番号同士で入れ替わる．
      do i=1, g_ipn-1  !i番目とi+1番目の隣り合う数値の節点座標の値を比較して，i番目の方が値が大きければ番号を入れ替える
        if(g_node(g_element(a,th_ge(i)), 3) < g_node(g_element(a,th_ge(i+1)), 3))then               !global節点座標のx成分の大小を比較
          temp = 0
          temp = th_ge(i)
          th_ge(i) = th_ge(i+1)
          th_ge(i+1) = temp
        end if
      end do
    end do


   !2.th_ge(1)～th_ge(4),th_ge(5)～th_ge(8)それぞれでy座標の大きいものから順に並び替え
    !th_ge(1)～th_ge(4)内でz座標の大小で並び替え
    do k=1, (g_ipn-4)-1  !昇順が完全に完了させるためにglobal_ipn-1回繰り返す．1回の操作で隣り合う番号同士で入れ替わる．
      do i=1, (g_ipn-4)-1  !i番目とi+1番目の隣り合う数値の節点座標の値を比較して，i番目の方が値が大きければ番号を入れ替える
        if(g_node(g_element(a,th_ge(i)), 2) < g_node(g_element(a,th_ge(i+1)), 2))then               !global節点座標のy成分の大小を比較
          temp = 0
          temp = th_ge(i)
          th_ge(i) = th_ge(i+1)
          th_ge(i+1) = temp
        end if
      end do
    end do

   !th_ge(5)～th_ge(8)内でz座標の大小で並び替え
    do k=1, (g_ipn-4)-1  !昇順が完全に完了させるためにglobal_ipn-1回繰り返す．1回の操作で隣り合う番号同士で入れ替わる．
      do i=5, g_ipn-1  !i番目とi+1番目の隣り合う数値の節点座標の値を比較して，i番目の方が値が大きければ番号を入れ替える
        if(g_node(g_element(a,th_ge(i)), 2) < g_node(g_element(a,th_ge(i+1)), 2))then               !global節点座標のy成分の大小を比較
          temp = 0
          temp = th_ge(i)
          th_ge(i) = th_ge(i+1)
          th_ge(i+1) = temp
        end if
      end do
    end do
    

      !3.y座標の大きい4節点と小さい4節点それぞれの中央点の座標を計算し，要素の傾きを求める
      node_e1 = 0.0d0
      node_e2 = 0.0d0
      do k=1, 2 !x,y座標のみ計算
        node_e1(k) = (g_node(g_element(a,th_ge(1)),k) &
                   + g_node(g_element(a,th_ge(2)),k) &
                   + g_node(g_element(a,th_ge(5)),k) &
                   + g_node(g_element(a,th_ge(6)),k)) / 4
        node_e2(k) = (g_node(g_element(a,th_ge(3)),k) &
                   + g_node(g_element(a,th_ge(4)),k) &
                   + g_node(g_element(a,th_ge(7)),k) &
                   + g_node(g_element(a,th_ge(8)),k)) / 4
      end do


 
      !-----------------------------------------------------------------各要素の繊維方向計算-----------------------------------------------------------------------------------------------------------------
      !z軸中心の回転角を計算
      th(a) = atan((node_e1(1) - node_e2(1)) / (node_e1(2) - node_e2(2))) !y軸基準の傾き角
      
      

      !今回のモデル（全部の節点がプラスの値で，z方向に厚みを持たせたダブテールモデル）の場合の補正
      temp_th(a)=180*(-1 * th(a))/pi             !今回は繊維方向の傾き角がxy平面上で傾きがマイナスとなるので"-1"を掛けている
      if(abs(temp_th(a)) <= 0.01 )then
        temp_th(a) = 0.0d0
      end if
print *, 'calculated 4'


      !傾きが大きいかつdx-dyのアスペクト比が小さい要素のはじき出し
      if(temp_th(a) < -10.0d0 .or. temp_th(a) > 50.0d0) then
        temp = th_ge(2)
        th_ge(2) = th_ge(3)
        th_ge(3) = temp
        temp = th_ge(6)
        th_ge(6) = th_ge(7)
        th_ge(7) = temp
      end if
print *, 'calculated 5'



  !2.y座標の大きい4節点と小さい4節点それぞれの中央点の座標を計算し，要素の傾きを求める
      node_e1 = 0.0d0
      node_e2 = 0.0d0
      do k=1, 2 !x,y座標のみ計算
       node_e1(k) = (g_node(g_element(a,th_ge(1)),k) &
                   + g_node(g_element(a,th_ge(2)),k) &
                   + g_node(g_element(a,th_ge(5)),k) &
                   + g_node(g_element(a,th_ge(6)),k)) / 4
        node_e2(k) = (g_node(g_element(a,th_ge(3)),k) &
                   + g_node(g_element(a,th_ge(4)),k) &
                   + g_node(g_element(a,th_ge(7)),k) &
                   + g_node(g_element(a,th_ge(8)),k)) / 4
      end do
print *, 'calculated 6'

      !z軸中心の回転角を再計算
      th(a) = atan((node_e1(1) - node_e2(1)) / (node_e1(2) - node_e2(2))) !y軸基準の傾き角

      !今回のモデル（全部の節点がプラスの値で，x方向に厚みを持たせたダブテールモデル）の場合の補正
      temp_th(a)=180*(-1 * th(a))/pi             !今回は繊維方向の傾き角がyz平面上で傾きがマイナスとなるので"-1"を掛けている
      if(abs(temp_th(a)) <= 0.01 )then
        temp_th(a) = 0.0d0
      end if

      !th(a)の値に戻す
      th(a) = pi*temp_th(a)/180    
      print *, 'calculated 6'
      
      end do !aの終わり
      
     !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      !角度計算はここまで．以下sorting_element用

   
    
  !繊維の傾き角度の出力
 
    open(10,file='output_th.dat')
      do i=1,g_econst
          write(10,'(I6,7F15.5)') i,temp_th(i)
      end do
    close(10)
      
    
    open(11,file='output_check_th.dat')
      do i=1, g_econst
        write(11,'(I6, 7F15.5)') i,check_th(i,1),check_th(i,2) 
      end do
    close(11)
  
   
  
  
  
    print *, 'make theta'
    
  end subroutine
    