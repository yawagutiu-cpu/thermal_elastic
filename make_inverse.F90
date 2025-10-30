subroutine make_inverse(A, invA, m)        !LAPACKを用いて逆行列を求める関数

use parameter1
implicit none

integer n, i, j
integer, intent(in) :: m       !intent: inなら入力用引数，outなら出力用引数，inoutならどっちも を指定
double precision, intent(in) :: A(m, m)
double precision, intent(inout) :: invA(m, m)

double precision :: keepA(m, m)
double precision :: invA_temp(m, m)

!各種パラメータ初期化
ipiv      = 0
info      = 0
keepA     = 0.0d0
invA_temp = 0.0d0
invA      = 0.0d0

!keepAにAの行列を記憶させておいて，計算後元のAが消えないようにする
keepA = A

!invAを単位行列にする
do n=1, m
    invA_temp(n,n) = 1.0d0
end do

!LAPACKを実行
call dgetrf(m, m, keepA, m, ipiv, info)
!           (m, n, keepA, lda, ipiv, info)
call dgetrs(trans, m, m, keepA, m, ipiv, invA_temp, m, info)
!           (trans, n, nrhs, keepA, lda, ipiv, X, ldb, info)

invA = invA_temp


end subroutine