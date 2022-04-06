module LibMdoodz

using mdoodz_jll
export mdoodz_jll

using CEnum

struct BC
    type::Ptr{Cchar}
    val::Ptr{Cdouble}
end

struct BCT
    type::Ptr{Cchar}
    typW::Ptr{Cchar}
    typE::Ptr{Cchar}
    typS::Ptr{Cchar}
    typN::Ptr{Cchar}
    val::Ptr{Cdouble}
    valW::Ptr{Cdouble}
    valE::Ptr{Cdouble}
    valS::Ptr{Cdouble}
    valN::Ptr{Cdouble}
end

struct grid
    Nx::Cint
    Nz::Cint
    NN::Cint
    NC::Cint
    dx::Cdouble
    dz::Cdouble
    roger_x::Ptr{Cdouble}
    roger_z::Ptr{Cdouble}
    div_u::Ptr{Cdouble}
    div_u_s::Ptr{Cdouble}
    div_u_el::Ptr{Cdouble}
    div_u_pl::Ptr{Cdouble}
    div_u_r::Ptr{Cdouble}
    u_in::Ptr{Cdouble}
    v_in::Ptr{Cdouble}
    p_in::Ptr{Cdouble}
    p_corr::Ptr{Cdouble}
    sxxd::Ptr{Cdouble}
    szzd::Ptr{Cdouble}
    sxz::Ptr{Cdouble}
    exxd::Ptr{Cdouble}
    ezzd::Ptr{Cdouble}
    exz::Ptr{Cdouble}
    VE_s::Ptr{Cdouble}
    VE_n::Ptr{Cdouble}
    sxxd0::Ptr{Cdouble}
    szzd0::Ptr{Cdouble}
    sxz0::Ptr{Cdouble}
    mu_s::Ptr{Cdouble}
    mu_n::Ptr{Cdouble}
    u_adv::Ptr{Cdouble}
    v_adv::Ptr{Cdouble}
    eta_phys_n::Ptr{Cdouble}
    kx::Ptr{Cdouble}
    kz::Ptr{Cdouble}
    Cv::Ptr{Cdouble}
    Qr::Ptr{Cdouble}
    eta_phys_s::Ptr{Cdouble}
    u_start::Ptr{Cdouble}
    v_start::Ptr{Cdouble}
    p_start::Ptr{Cdouble}
    divth0_n::Ptr{Cdouble}
    T0_n::Ptr{Cdouble}
    iter_smooth::Ptr{Cint}
    nb_part_cell::Ptr{Cint}
    nb_part_vert::Ptr{Cint}
    BCu::BC
    BCv::BC
    BCp::BC
    BCp_exp::BC
    BCt::BCT
    BCg::BCT
    BCc::BCT
    xg_coord::Ptr{Cdouble}
    zg_coord::Ptr{Cdouble}
    xc_coord::Ptr{Cdouble}
    zc_coord::Ptr{Cdouble}
    xvz_coord::Ptr{Cdouble}
    zvx_coord::Ptr{Cdouble}
    xg_coord0::Ptr{Cdouble}
    zg_coord0::Ptr{Cdouble}
    xg_coord_ext::Ptr{Cdouble}
    zg_coord_ext::Ptr{Cdouble}
    eta_s::Ptr{Cdouble}
    eta_n::Ptr{Cdouble}
    rho_s::Ptr{Cdouble}
    rho_n::Ptr{Cdouble}
    X_s::Ptr{Cdouble}
    X_n::Ptr{Cdouble}
    X0_s::Ptr{Cdouble}
    X0_n::Ptr{Cdouble}
    p0_n::Ptr{Cdouble}
    p0_s::Ptr{Cdouble}
    OverS_n::Ptr{Cdouble}
    strain_n::Ptr{Cdouble}
    strain_s::Ptr{Cdouble}
    u::Ptr{Cdouble}
    v::Ptr{Cdouble}
    p::Ptr{Cdouble}
    ru::Ptr{Cdouble}
    rv::Ptr{Cdouble}
    rp::Ptr{Cdouble}
    rhs_u::Ptr{Cdouble}
    rhs_v::Ptr{Cdouble}
    rhs_p::Ptr{Cdouble}
    rhs_t::Ptr{Cdouble}
    gx::Ptr{Cdouble}
    gz::Ptr{Cdouble}
    p_scale::Cdouble
    alp::Ptr{Cdouble}
    bet_n::Ptr{Cdouble}
    bet_s::Ptr{Cdouble}
    p_lith::Ptr{Cdouble}
    p_lith0::Ptr{Cdouble}
    dp::Ptr{Cdouble}
    Qrho::Ptr{Cdouble}
    VxVz::Ptr{Cdouble}
    VzVx::Ptr{Cdouble}
    P2N::Ptr{Cint}
    P2C::Ptr{Cint}
    kvx::Ptr{Cint}
    lvx::Ptr{Cint}
    kvz::Ptr{Cint}
    lvz::Ptr{Cint}
    kp::Ptr{Cint}
    lp::Ptr{Cint}
    kn::Ptr{Cint}
    ln::Ptr{Cint}
    phase_perc_n::Ptr{Ptr{Cdouble}}
    phase_perc_s::Ptr{Ptr{Cdouble}}
    phase_eta_n::Ptr{Ptr{Cdouble}}
    phase_eta_s::Ptr{Ptr{Cdouble}}
    sxxd0_s::Ptr{Cdouble}
    szzd0_s::Ptr{Cdouble}
    sxz0_n::Ptr{Cdouble}
    exxd_s::Ptr{Cdouble}
    ezzd_s::Ptr{Cdouble}
    exz_n::Ptr{Cdouble}
    sxz_n::Ptr{Cdouble}
    rho0_n::Ptr{Cdouble}
    Uthermal::Cdouble
    Uelastic::Cdouble
    Work::Cdouble
    Tii_mean::Cdouble
    Eii_mean::Cdouble
    T_mean::Cdouble
    P_mean::Cdouble
    Work_time::Ptr{Cdouble}
    Uelastic_time::Ptr{Cdouble}
    Uthermal_time::Ptr{Cdouble}
    Time_time::Ptr{Cdouble}
    Short_time::Ptr{Cdouble}
    P_mean_time::Ptr{Cdouble}
    T_mean_time::Ptr{Cdouble}
    Tii_mean_time::Ptr{Cdouble}
    Eii_mean_time::Ptr{Cdouble}
    T::Ptr{Cdouble}
    dT::Ptr{Cdouble}
    d_n::Ptr{Cdouble}
    d0_n::Ptr{Cdouble}
    phi_n::Ptr{Cdouble}
    phi0_n::Ptr{Cdouble}
    eII_el::Ptr{Cdouble}
    eII_pl::Ptr{Cdouble}
    eII_pl_s::Ptr{Cdouble}
    eII_pwl::Ptr{Cdouble}
    eII_exp::Ptr{Cdouble}
    eII_lin::Ptr{Cdouble}
    eII_gbs::Ptr{Cdouble}
    eII_cst::Ptr{Cdouble}
    eII_pwl_s::Ptr{Cdouble}
    exx_el::Ptr{Cdouble}
    ezz_el::Ptr{Cdouble}
    exz_el::Ptr{Cdouble}
    exx_diss::Ptr{Cdouble}
    ezz_diss::Ptr{Cdouble}
    exz_diss::Ptr{Cdouble}
    comp_cells::Ptr{Cint}
    D11_n::Ptr{Cdouble}
    D12_n::Ptr{Cdouble}
    D13_n::Ptr{Cdouble}
    D14_n::Ptr{Cdouble}
    D21_n::Ptr{Cdouble}
    D22_n::Ptr{Cdouble}
    D23_n::Ptr{Cdouble}
    D24_n::Ptr{Cdouble}
    D31_s::Ptr{Cdouble}
    D32_s::Ptr{Cdouble}
    D33_s::Ptr{Cdouble}
    D34_s::Ptr{Cdouble}
    detadexx_n::Ptr{Cdouble}
    detadezz_n::Ptr{Cdouble}
    detadgxz_n::Ptr{Cdouble}
    detadp_n::Ptr{Cdouble}
    ddivpdexx_n::Ptr{Cdouble}
    ddivpdezz_n::Ptr{Cdouble}
    ddivpdgxz_n::Ptr{Cdouble}
    ddivpdp_n::Ptr{Cdouble}
    detadexx_s::Ptr{Cdouble}
    detadezz_s::Ptr{Cdouble}
    detadgxz_s::Ptr{Cdouble}
    detadp_s::Ptr{Cdouble}
    drhodp_n::Ptr{Cdouble}
    phi0_s::Ptr{Cdouble}
    d0_s::Ptr{Cdouble}
    T_s::Ptr{Cdouble}
    P_s::Ptr{Cdouble}
    FS_AR_n::Ptr{Cdouble}
    FS_AR_s::Ptr{Cdouble}
    aniso_factor_n::Ptr{Cdouble}
    aniso_factor_s::Ptr{Cdouble}
    d1_n::Ptr{Cdouble}
    d2_n::Ptr{Cdouble}
    d1_s::Ptr{Cdouble}
    d2_s::Ptr{Cdouble}
    cell_min_z::Ptr{Cdouble}
    cell_max_z::Ptr{Cdouble}
    vert_min_z::Ptr{Cdouble}
    vert_max_z::Ptr{Cdouble}
    dil_n::Ptr{Cdouble}
    dil_s::Ptr{Cdouble}
    fric_n::Ptr{Cdouble}
    fric_s::Ptr{Cdouble}
    C_n::Ptr{Cdouble}
    C_s::Ptr{Cdouble}
    exz_n_el::Ptr{Cdouble}
    exz_n_diss::Ptr{Cdouble}
    exz_n_pl::Ptr{Cdouble}
    Wdiss::Ptr{Cdouble}
    Wel::Ptr{Cdouble}
    Wtot::Ptr{Cdouble}
    kc_x::Ptr{Cdouble}
    kc_z::Ptr{Cdouble}
    FreeSurfW_s::Ptr{Cdouble}
    FreeSurfW_n::Ptr{Cdouble}
    noise_n::Ptr{Cdouble}
    noise_s::Ptr{Cdouble}
end

@cenum ETA_AVG::UInt32 begin
    ARITHMETIC = 0
    HARMONIC = 1
    GEOMETRIC = 2
end

struct params
    xmin::Cdouble
    zmin::Cdouble
    xmax::Cdouble
    zmax::Cdouble
    time::Cdouble
    dx::Cdouble
    dz::Cdouble
    dt::Cdouble
    dt0::Cdouble
    dt_start::Cdouble
    dt_max::Cdouble
    L0::Cdouble
    dt_min::Cdouble
    xmin0::Cdouble
    zmin0::Cdouble
    xmax0::Cdouble
    zmax0::Cdouble
    gx::Cdouble
    gz::Cdouble
    Nx::Cint
    Nz::Cint
    Nt::Cint
    step::Cint
    nit::Cint
    Newton::Cint
    noisy::Cint
    eta_avg::ETA_AVG
    itp_stencil::Cint
    nexp_radial_basis::Cdouble
    ismechanical::Cint
    isperiodic_x::Cint
    isinertial::Cint
    iselastic::Cint
    isnonnewtonian::Cint
    isthermal::Cint
    ispureshear_ale::Cint
    free_surf::Cint
    write_markers::Cint
    write_debug::Cint
    free_surf_stab::Cdouble
    dt_constant::Cint
    RK::Cint
    line_search::Cint
    thermal_eq::Cint
    subgrid_diff::Cint
    adiab_heat::Cint
    shear_heat::Cint
    advection::Cint
    fstrain::Cint
    ConservInterp::Cint
    surf_processes::Cint
    cpc::Cint
    surf_remesh::Cint
    loc_iter::Cint
    therm_pert::Cint
    surf_ised1::Cint
    surf_ised2::Cint
    MantleID::Cint
    topografix::Cint
    Reseed::Cint
    SmoothSoftening::Cint
    EpsBG::Cdouble
    DivBG::Cdouble
    user0::Cdouble
    user1::Cdouble
    user2::Cdouble
    user3::Cdouble
    user4::Cdouble
    user5::Cdouble
    user6::Cdouble
    user7::Cdouble
    user8::Cdouble
    input_file::Ptr{Cchar}
    Nb_phases::Cint
    ncont::Cint
    Courant::Cdouble
    mineta::Cdouble
    maxeta::Cdouble
    initial_noise::Cint
    initial_part::Cint
    decoupled_solve::Cint
    lsolver::Cint
    diag_scaling::Cint
    pc_type::Cint
    penalty::Cdouble
    abs_tol_div::Cdouble
    rel_tol_div::Cdouble
    auto_penalty::Cdouble
    compressible::Cdouble
    rel_tol_KSP::Cdouble
    line_search_min::Cdouble
    safe_dt_div::Cdouble
    safe_mode::Cint
    nstagmax::Cint
    nT::Cint
    nE::Cint
    nd::Cint
    def_maps::Cint
    Pn::Cdouble
    Tmin::Cdouble
    Tmax::Cdouble
    Emin::Cdouble
    Emax::Cdouble
    dmin::Cdouble
    dmax::Cdouble
    PrBG::Cdouble
    TBG::Cdouble
    surf_diff::Cdouble
    surf_sedirate::Cdouble
    surf_baselev::Cdouble
    surf_Winc::Cdouble
    surf_Vinc::Cdouble
    therm_pert_x0::Cdouble
    therm_pert_z0::Cdouble
    therm_pert_dT::Cdouble
    therm_pert_rad::Cdouble
    cooling_time::Cdouble
    force_act_vol_ast::Cint
    act_vol_dis_ast::Cdouble
    act_vol_dif_ast::Cdouble
    isPD::Cint
    num_PD::Cint
    PDMnT::Ptr{Cint}
    PDMnP::Ptr{Cint}
    PD1DnP::Ptr{Cint}
    PDMrho::Ptr{Ptr{Cdouble}}
    PDMTmin::Ptr{Cdouble}
    PDMTmax::Ptr{Cdouble}
    PDMPmin::Ptr{Cdouble}
    PDMPmax::Ptr{Cdouble}
    PD1Drho::Ptr{Ptr{Cdouble}}
    PD1Dmin::Ptr{Cdouble}
    PD1Dmax::Ptr{Cdouble}
    kin_nP::Cint
    kin_nT::Cint
    kin_dG::Ptr{Cdouble}
    kin_Tmin::Cdouble
    kin_Tmax::Cdouble
    kin_Pmin::Cdouble
    kin_Pmax::Cdouble
    rec_T_P_x_z::Cint
    delete_breakpoints::Cint
    GNUplot_residuals::Cint
    BC_setup_type::Cint
    shear_style::Cint
    polar::Cint
    StressRotation::Cint
    StressUpdate::Cint
    DirectNeighbour::Cint
    diffuse_X::Cint
    diffuse_avg::Cint
    diffusion_length::Cdouble
    ProgReac::Cint
    NoReturn::Cint
    VolChangeReac::Cint
    Plith_trick::Cint
    UnsplitDiffReac::Cint
    kinetics::Cint
    aniso::Cint
    aniso_fstrain::Cint
    oop::Cint
    noise_bg::Cint
    eqn_state::Cint
    residual_form::Cint
end

struct scale
    eta::Cdouble
    L::Cdouble
    V::Cdouble
    T::Cdouble
    t::Cdouble
    a::Cdouble
    E::Cdouble
    S::Cdouble
    m::Cdouble
    rho::Cdouble
    F::Cdouble
    J::Cdouble
    W::Cdouble
    Cv::Cdouble
    rhoE::Cdouble
    k::Cdouble
end

struct markers
    Nx_part::Cint
    Nz_part::Cint
    Nb_part::Cint
    Nb_part_max::Cint
    min_part_cell::Cint
    Nb_part_ini::Cint
    x::Ptr{Cdouble}
    z::Ptr{Cdouble}
    Vx::Ptr{Cdouble}
    Vz::Ptr{Cdouble}
    P::Ptr{Cdouble}
    sxxd::Ptr{Cdouble}
    szzd::Ptr{Cdouble}
    sxz::Ptr{Cdouble}
    progress::Ptr{Cdouble}
    T::Ptr{Cdouble}
    d::Ptr{Cdouble}
    phi::Ptr{Cdouble}
    X::Ptr{Cdouble}
    syy::Ptr{Cdouble}
    dsyy::Ptr{Cdouble}
    strain::Ptr{Cdouble}
    strain_el::Ptr{Cdouble}
    strain_pl::Ptr{Cdouble}
    strain_pwl::Ptr{Cdouble}
    strain_exp::Ptr{Cdouble}
    strain_lin::Ptr{Cdouble}
    strain_gbs::Ptr{Cdouble}
    phase::Ptr{Cint}
    generation::Ptr{Cint}
    dual::Ptr{Cint}
    intag::Ptr{Cint}
    Fxx::Ptr{Cdouble}
    Fxz::Ptr{Cdouble}
    Fzx::Ptr{Cdouble}
    Fzz::Ptr{Cdouble}
    nx::Ptr{Cdouble}
    nz::Ptr{Cdouble}
    T0::Ptr{Cdouble}
    P0::Ptr{Cdouble}
    x0::Ptr{Cdouble}
    z0::Ptr{Cdouble}
    Tmax::Ptr{Cdouble}
    Pmax::Ptr{Cdouble}
    divth::Ptr{Cdouble}
    dsxxd::Ptr{Cdouble}
    dszzd::Ptr{Cdouble}
    dsxz::Ptr{Cdouble}
    noise::Ptr{Cdouble}
    rho::Ptr{Cdouble}
end

struct mat_prop
    Nb_phases::Cint
    R::Cdouble
    eta0::NTuple{20, Cdouble}
    rho::NTuple{20, Cdouble}
    mu::NTuple{20, Cdouble}
    Cv::NTuple{20, Cdouble}
    k::NTuple{20, Cdouble}
    Qr::NTuple{20, Cdouble}
    C::NTuple{20, Cdouble}
    phi::NTuple{20, Cdouble}
    psi::NTuple{20, Cdouble}
    Slim::NTuple{20, Cdouble}
    n::NTuple{20, Cdouble}
    A::NTuple{20, Cdouble}
    Ea::NTuple{20, Cdouble}
    Va::NTuple{20, Cdouble}
    alp::NTuple{20, Cdouble}
    bet::NTuple{20, Cdouble}
    Qm::NTuple{20, Cdouble}
    T0::NTuple{20, Cdouble}
    P0::NTuple{20, Cdouble}
    drho::NTuple{20, Cdouble}
    k_eff::NTuple{20, Cdouble}
    tpwl::NTuple{20, Cdouble}
    Qpwl::NTuple{20, Cdouble}
    Vpwl::NTuple{20, Cdouble}
    npwl::NTuple{20, Cdouble}
    mpwl::NTuple{20, Cdouble}
    Apwl::NTuple{20, Cdouble}
    apwl::NTuple{20, Cdouble}
    fpwl::NTuple{20, Cdouble}
    rpwl::NTuple{20, Cdouble}
    Fpwl::NTuple{20, Cdouble}
    pref_pwl::NTuple{20, Cdouble}
    texp::NTuple{20, Cdouble}
    Qexp::NTuple{20, Cdouble}
    Vexp::NTuple{20, Cdouble}
    Sexp::NTuple{20, Cdouble}
    Eexp::NTuple{20, Cdouble}
    Gexp::NTuple{20, Cdouble}
    aexp::NTuple{20, Cdouble}
    fexp::NTuple{20, Cdouble}
    rexp::NTuple{20, Cdouble}
    qexp::NTuple{20, Cdouble}
    nexp::NTuple{20, Cdouble}
    tlin::NTuple{20, Cdouble}
    Qlin::NTuple{20, Cdouble}
    Vlin::NTuple{20, Cdouble}
    nlin::NTuple{20, Cdouble}
    mlin::NTuple{20, Cdouble}
    Alin::NTuple{20, Cdouble}
    alin::NTuple{20, Cdouble}
    flin::NTuple{20, Cdouble}
    rlin::NTuple{20, Cdouble}
    Flin::NTuple{20, Cdouble}
    tgbs::NTuple{20, Cdouble}
    Qgbs::NTuple{20, Cdouble}
    Vgbs::NTuple{20, Cdouble}
    ngbs::NTuple{20, Cdouble}
    mgbs::NTuple{20, Cdouble}
    Agbs::NTuple{20, Cdouble}
    agbs::NTuple{20, Cdouble}
    fgbs::NTuple{20, Cdouble}
    rgbs::NTuple{20, Cdouble}
    Fgbs::NTuple{20, Cdouble}
    ppzm::NTuple{20, Cdouble}
    Kpzm::NTuple{20, Cdouble}
    Qpzm::NTuple{20, Cdouble}
    Vpzm::NTuple{20, Cdouble}
    Gpzm::NTuple{20, Cdouble}
    cpzm::NTuple{20, Cdouble}
    Lpzm::NTuple{20, Cdouble}
    gs_ref::NTuple{20, Cdouble}
    Skin::NTuple{20, Cdouble}
    kkin::NTuple{20, Cdouble}
    Qkin::NTuple{20, Cdouble}
    kin::NTuple{20, Cint}
    gs::NTuple{20, Cint}
    cstv::NTuple{20, Cint}
    pwlv::NTuple{20, Cint}
    linv::NTuple{20, Cint}
    expv::NTuple{20, Cint}
    gbsv::NTuple{20, Cint}
    phase_diagram::NTuple{20, Cint}
    density_model::NTuple{20, Cint}
    C_end::NTuple{20, Cdouble}
    phi_end::NTuple{20, Cdouble}
    psi_end::NTuple{20, Cdouble}
    pls_start::NTuple{20, Cdouble}
    pls_end::NTuple{20, Cdouble}
    eta_vp::NTuple{20, Cdouble}
    n_vp::NTuple{20, Cdouble}
    plast::NTuple{20, Cint}
    phi_soft::NTuple{20, Cint}
    psi_soft::NTuple{20, Cint}
    coh_soft::NTuple{20, Cint}
    is_tensile::NTuple{20, Cint}
    Pr::NTuple{20, Cdouble}
    tau_kin::NTuple{20, Cdouble}
    dPr::NTuple{20, Cdouble}
    k_chem::NTuple{20, Cdouble}
    reac_soft::NTuple{20, Cint}
    reac_phase::NTuple{20, Cint}
    phase_mix::NTuple{20, Cint}
    phase_two::NTuple{20, Cint}
    aniso_factor::NTuple{20, Cdouble}
    aniso_angle::NTuple{20, Cdouble}
end

struct surface
    a::Ptr{Cdouble}
    b::Ptr{Cdouble}
    height::Ptr{Cdouble}
    vx::Ptr{Cdouble}
    vz::Ptr{Cdouble}
    a0::Ptr{Cdouble}
    b0::Ptr{Cdouble}
    height0::Ptr{Cdouble}
    VertInd::Ptr{Cint}
end

function GetSetupFileName(nargs, args)
    ccall((:GetSetupFileName, libmdoodz), Ptr{Cchar}, (Cint, Ptr{Ptr{Cchar}}), nargs, args)
end

function MinMaxArray(array, scale_, size, text)
    ccall((:MinMaxArray, libmdoodz), Cvoid, (Ptr{Cdouble}, Cdouble, Cint, Ptr{Cchar}), array, scale_, size, text)
end

# typedef void ( * BuildInitialTopography ) ( markers * topo_chain , params model , scale scaling )
const BuildInitialTopography = Ptr{Cvoid}

# typedef void ( * SetParticles ) ( markers * particles , scale scaling , params model , mat_prop * materials )
const SetParticles = Ptr{Cvoid}

# typedef void ( * SetBCs ) ( grid * mesh , params * model , scale scaling , markers * particles , mat_prop * materials , surface * topo )
const SetBCs = Ptr{Cvoid}

function RunMDOODZ(inputFileName, arg2, arg3, arg4)
    ccall((:RunMDOODZ, libmdoodz), Cint, (Ptr{Cchar}, BuildInitialTopography, SetParticles, SetBCs), inputFileName, arg2, arg3, arg4)
end

const zeroC = 273.15

const Rg = 8.31451

const PI = 3.14159265359

const Rad_Earth = 6370000

# exports
const PREFIXES = ["CX", "clang_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
