subroutine B_coordinates
!ガウス点座標、ξηζ座標を設定 何も変わらず

use parameter1
implicit none

!配列サイズを指定
allocate (gauss(8,3))
allocate (polar(8,3))

!配列初期化
gauss = 0.0d0
polar = 0.0d0

!ガウス点座標
gauss(1, 1:3) = (/-ga, -ga, -ga/)
gauss(2, 1:3) = (/ ga, -ga, -ga/)
gauss(3, 1:3) = (/ ga,  ga, -ga/)
gauss(4, 1:3) = (/-ga,  ga, -ga/)
gauss(5, 1:3) = (/-ga, -ga,  ga/)
gauss(6, 1:3) = (/ ga, -ga,  ga/)
gauss(7, 1:3) = (/ ga,  ga,  ga/)
gauss(8, 1:3) = (/-ga,  ga,  ga/)

!直交座標のξ-η-ζ座標系における座標
polar(1, 1:3) = (/-1.0d0, -1.0d0, -1.0d0/)
polar(2, 1:3) = (/ 1.0d0, -1.0d0, -1.0d0/)
polar(3, 1:3) = (/ 1.0d0,  1.0d0, -1.0d0/)
polar(4, 1:3) = (/-1.0d0,  1.0d0, -1.0d0/)
polar(5, 1:3) = (/-1.0d0, -1.0d0,  1.0d0/)
polar(6, 1:3) = (/ 1.0d0, -1.0d0,  1.0d0/)
polar(7, 1:3) = (/ 1.0d0,  1.0d0,  1.0d0/)
polar(8, 1:3) = (/-1.0d0,  1.0d0,  1.0d0/)

end subroutine