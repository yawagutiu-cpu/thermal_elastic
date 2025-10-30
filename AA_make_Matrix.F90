subroutine AA_make_Matrix
!配列サイズを指定(inputで必要なもののみ)   とりあえずalphaを追加
use parameter1
implicit none


allocate(deg_f(4))
allocate(e_rate(5))
allocate(deg_flag(4))

allocate(ipiv(nconst*3))

allocate(node(nconst,3))
allocate(element(econst,ipn))
allocate(g_node(g_nconst,g_ipn))
allocate(g_element(g_econst,g_ipn))
allocate(pid(g_econst))
allocate(ply_deg(g_ply))



allocate(young(kind))
allocate(poisson(kind))
!熱弾性
allocate(alpha(kind+3,6))

allocate(material(econst))
allocate(fiber_dir(g_econst))

allocate(option_file(40))

end subroutine
