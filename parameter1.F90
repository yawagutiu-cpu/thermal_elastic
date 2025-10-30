module parameter1
!みかん　とりあえず1要素ツースケール熱弾性をコピペ、ダブテールに必要なものを後々追加していく
implicit none

!基本情報(basic_info)にて入力---------------------------------------------------------------
integer nconst         !節点数
integer econst         !要素数
integer ipn            !積分点数(形状関数の数)
integer dim            !次元数(y座標系)
integer mode           !解析種類, 方法
double precision vol  !ユニットセル体積
double precision ga   !ガウス点　√1/3, 重みwi=1
double precision pi   !円周率
double precision eps  !擾乱倍率
integer kind           !等方性材料の数

integer g_dim !グローバルモデルの次元(x座標系)
integer g_econst !グローバルモデルの要素数
integer g_nconst !グローバルモデルの節点数
integer g_ply
double precision g_vol !グローバルモデルの体積
integer g_ipn !グローバルモデルの積分点数
integer g_boundary !グローバルモデルの境界条件



integer C_penalty_flag !C対称性境界条件の設定スイッチ
integer Y_penalty_flag !Y対称性境界条件のスイッチ

integer g_Cpenalty_flag
integer g_Ypenalty_flag

integer trigger
integer g_mode
double precision dt !時間刻みのステップ
double precision l_mag !変位の倍率
double precision g_mag !グローバル変位の倍率

!以下、make_Matrixにてサイズを指定
integer,allocatable::element(:,:) !要素を構成する節点番号
double precision,allocatable::node(:,:) !節点座標
integer,allocatable::material(:) !要素ごとに材料の種類を指定

integer,allocatable::g_element(:,:) 
double precision,allocatable::g_node(:,:)
double precision,allocatable::fiber_dir(:) !繊維方向
double precision,allocatable::pid(:)  !その要素が何層目なのか
double precision,allocatable::ply_deg(:)

!負荷条件
integer dir !構造物の負荷方向
double precision g_end !構造物の荷重目標値
double precision g_giv !現在付加されている荷重
double precision g_load_rate !構造物の時間依存変形の荷重速度
double precision Umat_speed           !構造物の時間依存変形の変位速度[s^-1]
double precision l_temp

!材料定数
double precision,allocatable :: young(:)        !ヤング率(要素ごとに設定)
double precision,allocatable :: poisson(:)      !ポアソン比(要素ごとに設定)
double precision ELL                              !炭素繊維長手方向ヤング率
double precision ETT                              !炭素繊維径方向ヤング率
double precision VLL                              !炭素繊維ポアソン比ν31
double precision VTT                              !炭素繊維ポアソン比ν12
double precision GLT

!熱弾性
double precision deltaT   !温度変化
double precision alpha11  !異方性材料の繊維と直交方向線膨張係数
double precision alpha33 !異方性材料の繊維方向線膨張係数
double precision,allocatable:: alpha(:,:)  !線膨張係数
double precision,allocatable:: alpha_dummy(:)

double precision,allocatable :: Imat(:,:)    !単位行列

!FEM, 均質化で使用するパラメータ------------------------------------------------------------------------------

!LAPACK
integer,allocatable :: ipiv(:)
integer,allocatable :: ipiv_K(:)            !KmatをLU分解するのに必要
integer,allocatable :: ipiv_K_g(:)          !g_KmatをLU分解する用

character :: trans = "N"
integer info

!B_coordinates
double precision,allocatable::gauss(:,:) !Gauss点座標
double precision,allocatable::polar(:,:) !!ξηζ座標

!Local_Boundary_Condition ミクロな境界条件
double precision lambda   !ペナルティ数
integer pair               !ペアリング数
integer,allocatable :: pn(:,:)
integer,allocatable :: pn_C(:,:)
integer,allocatable :: pn_Y(:,:)
integer C1,C2,C3,C4

!--------------------------micro----------------------------------
!Dmat
double precision,allocatable::Dmat(:,:,:)
double precision,allocatable::keep_Dmat(:,:)
double precision,allocatable::keep_invDmat(:,:)
double precision,allocatable::invDmat(:,:,:)
double precision,allocatable::Smat(:,:)
double precision,allocatable::invSmat(:,:,:)
double precision,allocatable::invSmat_p(:,:)

!Bmat
double precision,allocatable::Bmat(:,:,:,:)
double precision,allocatable::Hmat(:,:,:)
double precision,allocatable::Jmat(:,:,:,:)
double precision,allocatable::invJmat(:,:,:,:)
double precision,allocatable::Jmat_p(:,:)
double precision,allocatable::invJmat_p(:,:)
double precision,allocatable::det_J(:,:)
double precision,allocatable::dNdxyz(:,:,:,:)

!Kmat
double precision,allocatable::tBDmat(:,:,:,:)
double precision,allocatable::tBDBmat(:,:,:)
double precision,allocatable::eKmat(:,:,:)
double precision,allocatable::Kmat(:,:)
double precision,allocatable::Kmat_LU(:,:)

!Fmat
double precision,allocatable::eFmat(:,:,:)
double precision,allocatable::Fmat(:,:)

!Xmat
double precision,allocatable::Xmat_temp(:,:)
double precision,allocatable::eXmat(:,:,:)
double precision,allocatable::Xmat(:,:)
double precision,allocatable::BXmat(:,:,:,:)

!Amat
double precision,allocatable::eAmat(:,:,:,:)
double precision,allocatable::Amat_temp(:,:)
double precision,allocatable::invAmat_temp(:,:)
!--------------------------熱弾性-----------------------------------------

!追加コンテンツ
double precision,allocatable :: l_S1(:,:,:)
double precision,allocatable :: l_S1e(:,:)
double precision,allocatable :: g_Se1(:,:)
double precision,allocatable :: g_Se2(:,:)
double precision,allocatable :: g_S1(:,:,:)


!make_Mmat
double precision,allocatable:: Mmat(:)
double precision,allocatable:: e_Mmat(:,:)

!make_Psi
double precision,allocatable:: Psi_vec(:)

!make_Gamma
double precision,allocatable:: e_Psi(:,:)
double precision,allocatable:: DB(:,:)
double precision,allocatable:: Dalpha(:,:,:)
double precision,allocatable:: DBPsi(:,:,:)
double precision,allocatable:: n_gamma(:,:,:) !積分点でのGamma
double precision,allocatable:: gamma(:)
double precision,allocatable:: macro_alpha(:)


!rotation_Gamma
double precision, allocatable :: rot_Gamma_ini(:,:)
double precision, allocatable :: Gamma_mat(:,:)
double precision, allocatable :: rot_Gamma_mat(:,:)
double precision, allocatable :: R_ini(:,:)
double precision, allocatable :: Rtensor(:,:,:)
double precision,allocatable :: g_rot_Gamma(:,:)
double precision,allocatable:: Qtemp(:,:,:)
double precision,allocatable:: Qrot_Gamma_ini(:,:,:)

!rotation_alpha
double precision,allocatable :: alpha_rot(:,:,:)

!make_BTgamma
!double precision,allocatable :: e_BTgamma(:,:,:)
double precision,allocatable :: e_BTgamma(:,:)
double precision,allocatable :: BTgamma(:)


!rotation_g_E
double precision,allocatable :: g_E_mat(:,:)
double precision,allocatable :: rot_g_E_ini(:,:)
double precision,allocatable :: Qrot_g_E_ini(:,:,:)
double precision,allocatable :: temp2(:,:,:)
double precision,allocatable :: rot_g_E(:,:,:)

!rot_g_E_th
double precision,allocatable :: rot_g_E_th(:,:,:)
double precision,allocatable :: rot_g_Ee_th(:,:)

!rot_g_Se_th
double precision,allocatable :: g_S_mat(:,:)
double precision,allocatable :: rot_g_S_ini(:,:)
double precision,allocatable :: rot_g_S_th(:,:,:)
double precision,allocatable :: rot_g_Se_th(:,:)



!--------------------------macro----------------------------------
!global_Bmat
double precision, allocatable ::g_Bmat(:,:,:,:)
double precision, allocatable ::g_Hmat(:,:,:)
double precision, allocatable ::g_Jmat(:,:,:,:)
double precision, allocatable ::g_invJmat(:,:,:,:)
double precision, allocatable ::g_det_J(:,:)
double precision, allocatable ::g_dNdxyz(:,:,:,:)

!global_Kmat
double precision, allocatable ::g_tBDmat(:,:,:,:)
double precision, allocatable ::g_eKmat(:,:,:)
double precision, allocatable ::g_Kmat(:,:)
double precision, allocatable ::g_Kmat_LU(:,:)
double precision, allocatable ::g_Kmat_dummy(:,:) !境界条件を適応する前のKmat

!グローバルモデルの境界条件に使用するパラメーター
double precision g_lambda
integer,allocatable ::g_pd(:,:)  !global_penalty_direction：固定する節点と固定方向を指定
!追加
double precision,allocatable :: x_spc(:)
  double precision,allocatable :: y_spc(:)
  double precision,allocatable :: z_spc(:)
  double precision,allocatable :: xyz_spc(:)
  double precision,allocatable :: z_up(:)
  integer spc_x_number
  integer spc_y_number
  integer spc_z_number
  integer z_up_number
integer spc_xyz_number

!時間依存で使用するパラメータ------------------------------------------------------------------------------
double precision ti
integer step !ステップ数
integer swich_end !最終ステップスイッチ
double precision disp_giv !与えるひずみ
double precision strain_giv !与える変位
double precision dist

!--------------------------micro----------------------------------
!Beta
double precision,allocatable::e_rate(:) !参照ひずみ速度
double precision,allocatable::Beta(:,:,:,:,:) !粘塑性関数
double precision,allocatable::hf(:,:) !硬化関数
double precision,allocatable::Ep_ac(:,:,:,:) !粘塑性ひずみ
double precision,allocatable::Ep_ac_e(:,:,:) !要素ごと粘塑性ひずみ
double precision,allocatable::S_d(:,:,:) !偏差応力
double precision,allocatable::S_e(:,:,:) !相当応力
double precision nonl        !粘塑性の非線形べき乗数
double precision hf_p1       !母材硬化関数の係数
double precision hf_p2       !母材硬化関数の係数
double precision hf_p3       !母材硬化関数の係数

!Gvec
double precision,allocatable::eGvec(:,:) !要素ごとの硬化関数
double precision,allocatable::Gvec(:) !硬化関数

!Pvec
double precision,allocatable::Gvec_temp(:) !硬化関数（方程式計算用）
double precision,allocatable::Pvec(:,:,:) !特性関数φ
double precision,allocatable::ePvec(:,:) !要素ごとの特性関数

!Rvec
double precision,allocatable::DBmat(:,:,:,:) !要素ごとのr-ベクトル
double precision,allocatable::eRvec(:,:,:) !要素ごとのr-ベクトル
double precision,allocatable::Rvec_s(:,:,:,:,:) !全積分点でのr-ベクトル
double precision,allocatable::Rvec_e(:,:,:,:)
double precision,allocatable::Rvec(:,:,:)
double precision,allocatable::tBRvec(:)
double precision,allocatable::etBRmat(:,:)

!---------------------FEM解析に使用するパラメータ----------------------------------------------------------------------------
!--------------------------Macro----------------------------------
!変位計算
double precision,allocatable::g_disp(:) !構造物の変位
double precision,allocatable::g_disp_speed(:) !構造物の変位速度
double precision,allocatable::g_disp_speed_tensor1(:,:) !変位速度を入れ替えてテンソル化したもの
double precision,allocatable::g_disp_speed_tensor2(:,:) !さらに入れ替え
double precision,allocatable::g_disp_speed_n(:) !仮想変位（非適合要素）
double precision,allocatable::g_node_after(:) !構造物の変形後の節点座標

!応力ひずみ計算

double precision,allocatable::disp_dummy1(:,:)
double precision,allocatable::disp_dummy2(:,:)

double precision,allocatable::strain_tensor(:,:) !ひずみテンソル
double precision,allocatable::g_S_speed(:,:,:) !応力速度
double precision,allocatable::g_E_speed(:,:,:) !ひずみ速度
double precision,allocatable::g_Se(:,:) !要素ごとの応力
double precision,allocatable::g_Ee(:,:) !要素ごとのひずみ
double precision,allocatable::g_S(:,:,:) !構造物の積分点におけるローカルなマクロ応力
double precision,allocatable::g_E(:,:,:) !構造物の積分点におけるローカルなマクロひずみ
double precision,allocatable::g_Beta(:,:,:)
double precision,allocatable::g_Beta_e(:,:)
double precision,allocatable::g_Ep_ac(:,:)
double precision,allocatable::g_Ep_ac_e(:)
double precision,allocatable::fvec(:) !節点に加わる力
double precision::total_load !要素全体に加わる力

!--------------------------micro----------------------------------
double precision,allocatable::l_S_speed(:,:,:) !ミクロ応力速度分布
double precision,allocatable::l_S(:,:,:,:,:) !積分点ごとのミクロ応力
double precision,allocatable::l_disp_sharp(:,:,:,:) !擾乱変位
double precision,allocatable::l_disp_sharp_speed(:) !擾乱変位速度
double precision,allocatable::l_node_after(:,:,:,:) !移動後の節点座標
double precision,allocatable::count_ns(:) !各節点のミクロ応力作成時のカウンター
double precision,allocatable::g_microstress_point(:,:) !各節点の応力
double precision,allocatable::g_micro_strain_point(:,:) !各節点のひずみ
double precision,allocatable::g_beta_point(:,:)
double precision,allocatable::beta_point(:,:,:,:)
double precision,allocatable::l_Se(:,:,:,:) !要素ごとのミクロ応力
double precision,allocatable::l_Ee(:,:) !要素ごとのミクロひずみ
double precision,allocatable::l_E(:,:,:,:) !ミクロひずみ（相当粘塑性ひずみ）
double precision,allocatable::Dmat_LU(:,:) !ミクロひずみを求めるときに使う
integer,allocatable::Ipiv_D(:) !DmatのLU分解に使う
double precision,allocatable::l_E_dummy(:)

!---------------------座標変換に使用するパラメータ----------------------------------------------------------------------------
!座標変換
integer rot_Amattion_plane !座標変換するトリガー
integer select !1=ひずみ制御 2=応力制御
double precision stress_x
double precision stress_y
double precision stress_shear
double precision stress
double precision strain
double precision,allocatable::Pmat(:,:) !マクロ応力回転行列
double precision,allocatable::invPmat(:,:) !逆行列
double precision,allocatable::Qmat(:,:) !マクロひずみ回転行列
double precision,allocatable::invQmat(:,:) !逆行列
double precision,allocatable::b_A(:,:) !座標変換に必要なコンプライアンス行列
double precision,allocatable::invb_A(:,:) !b_Aの逆行列
double precision,allocatable::QB(:,:)
double precision,allocatable::invPB(:,:) ![P][B]^-1
double precision :: Evm_rot
double precision :: Svm_rot !回転後のミーゼスひずみ、応力

!テンソル関連
double precision,allocatable::node_f(:,:) !nodeをx1方向に座標変換した後の節点座標
double precision,allocatable::node_dummy(:,:)
double precision,allocatable::rot_mat(:,:,:) !ミクロなユニットセルの回転行列
double precision,allocatable::rot_mat_rev(:,:,:) !逆回転して元に戻す回転行列
double precision,allocatable::node_rot(:,:,:) !node1-node3(deg_fごとに回転させたnode)
double precision,allocatable::Qmat_a(:,:,:) 
double precision,allocatable::Pmat_rev(:,:,:)
double precision,allocatable::rot_Amat(:,:,:,:)  ![P][B]^-1[Q]^-1
double precision,allocatable::tensor(:,:,:,:)
double precision,allocatable::rot_tensor(:,:,:,:)
double precision deg !座標回転角
double precision,allocatable::deg_f(:) !deg1-deg3
double precision,allocatable::deg_flag(:) !deg_flg1-deg_flag3


!rotationAmatのところ
double precision, allocatable :: g_Qmat(:, :)
double precision, allocatable :: g_Qmat_2(:, :)
double precision, allocatable :: g_Amat(:, :, :)
double precision, allocatable :: g_invAmat(:, :, :)
double precision,allocatable :: Amat_dummy(:,:)
double precision,allocatable:: Amat_inv_dummy(:,:)


!make_theata
double precision,allocatable:: th(:)   !要素の傾き（今回は３軸中心)   rad
integer,allocatable:: th_ge(:)
double precision,allocatable :: temp_th(:)   !度
double precision,allocatable:: check_th(:,:)
double precision,allocatable::th_single(:,:)

!---------------------NCE・EAS要素に使用するパラメータ----------------------------------------------------------------------------
!非適合要素
double precision,allocatable :: eKmat_NCE(:,:,:)            !非適合要素を考えた要素剛性マトリックス
double precision,allocatable :: eK_cc(:,:,:)                !適合要素分の要素剛性マトリックス
double precision,allocatable :: eK_cn(:,:,:)                !非適合要素分の要素剛性マトリックス
double precision,allocatable :: eK_nc(:,:,:)                !非適合要素分の要素剛性マトリックス
double precision,allocatable :: eK_nn(:,:,:)                !非適合要素分の要素剛性マトリックス
double precision,allocatable :: inveK_nn(:,:,:)             !非適合要素分の要素剛性マトリックスの逆行列
double precision,allocatable :: eKnn_keep(:,:)              !上記の逆行列を作るための一時的な非適合要素の要素剛性マトリックス
double precision,allocatable :: inveKnn_keep(:,:)           !上記の逆行列を作るための一時的な非適合要素の要素剛性マトリックスの逆行列
double precision,allocatable :: eKcn_inveKnn(:,:,:)         ![Keub]*[Kebb]^(-1)
double precision,allocatable :: eKcn_inveKnn_eKnc(:,:,:)    ![Keub]*[Kebb]^(-1)*[Kebu]
double precision,allocatable ::eFmat_NCE(:,:,:) 
double precision,allocatable :: eKcn_inveKnn_eFn(:,:,:)     ![Keub][Kebb]^(-1)*[Feb]
double precision,allocatable :: eF_c(:,:,:)                 !適合要素分の要素節点荷重
double precision,allocatable :: eF_n(:,:,:)                 !非適合要素分の要素節点荷重
double precision,allocatable :: eX_c(:,:,:)                 !適合要素分の特性関数
double precision,allocatable :: eX_n(:,:,:)                 !非適合要素分の特性関数
double precision,allocatable :: eKnc_eXc(:,:,:)             ![Kebu]*[Xe]

!EAS要素
double precision,allocatable :: eKmat_EAS(:,:,:)
double precision,allocatable :: eFmat_EAS(:,:,:)
double precision,allocatable :: H0mat(:,:)                  !Hマトリックス(拡張ひずみ成分)
double precision,allocatable :: J0mat(:,:,:)
double precision,allocatable :: detJ0(:)
double precision,allocatable :: Q0mat(:,:,:)                !拡張ひずみの座標変換行列
double precision,allocatable :: Q0mat_keep(:,:)
double precision,allocatable :: invQ0mat(:,:,:)
double precision,allocatable :: invQ0mat_keep(:,:)
double precision,allocatable :: Eep(:,:,:)                  !仮定されるひずみ
double precision,allocatable :: Gep(:,:,:,:)                !拡張ひずみ-変位マトリックス(ε~ = Gep * α)

!-------------------------------マクロ--------------------------------------------------------------------------------------------------------------------
double precision :: Evm                                  !ミーゼスひずみ!使ってない
double precision :: Svm                                  !ミーゼス応力

double precision,allocatable :: maS_speed_p(:)      !巨視的応力速度の部分ベクトル
double precision,allocatable :: maE_speed_p(:)      !巨視的ひずみ速度の部分ベクトル
double precision,allocatable :: maS_speed_rot(:)    !回転後の巨視的応力速度
double precision,allocatable :: maE_speed_rot(:)    !回転後の巨視的ひずみ速度
double precision,allocatable :: maS_speed_rot_p(:)  !回転後の巨視的応力速度ベクトル
double precision,allocatable :: maE_speed_rot_p(:)  !回転後の巨視的ひずみ速度ベクトル
double precision,allocatable :: maS_rot(:)          !回転後の巨視的応力
double precision,allocatable :: maE_rot(:)          !回転後の巨視的ひずみ
double precision,allocatable :: maE_tensor(:,:)     !巨視的ひずみを3*3テンソルに再定義したもの
double precision,allocatable :: BR(:)               ![A]^(-1)*{r}
double precision,allocatable :: BR_p(:)             ![A]^(-1)*{r}の部分ベクトル
double precision,allocatable :: QBR_p(:)            ![Q]*[B]*[R]

!-----------------データ入出力---------------------------------------------------------------------------------------------------------------------------
integer phase                                
character*50 fmt                                   !出力形式
character*50 name
character*200 filename
character*150 command
character*6 out_ele
integer file_number
integer cupling_number
integer step_number
integer mark
character*5 CWORKs                                 !ステップの文字変数
character*6 CWORKe                                 !要素の文字変数
character*2 CWORKi                                 !積分点の文字変数
integer,allocatable :: option_file(:)
integer,allocatable :: option_output_Ep(:,:)                !出力するマクロ関係の組み合わせ
double precision,allocatable :: option_output_step(:)       !出力するステップの情報
integer file_end

double precision :: x_temp, y_temp, z_temp
double precision :: x_middle, y_middle
integer, allocatable :: g_pn_swich(:, :)             !input_global_penaltyの1or0を格納

end module
