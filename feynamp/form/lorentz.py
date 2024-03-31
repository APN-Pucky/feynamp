from feynamp.form import *

polsum_phys = """
*id polsum(mu?,nu?) = -d_(mu,nu) + (pa(mu)*pb(nu)+pa(nu)*pb(mu))/pa.pb;
*id polsum(mu?,nu?) = -d_(mu,nu) + ax*((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);
id polsum(mu?,nu?) = -d_(mu,nu) + ((q(mu)*pb(nu)+q(nu)*pb(mu))/q.pb - q.q*pb(mu)*pb(nu)/q.pb/q.pb);
"""

gammas = """
repeat;
* identity
    id Gamma(Mua?,Spinb?,Spinc?) * GammaId(Spinc?,Spind?) = Gamma(Mua,Spinb,Spind);
    id Gamma(Mua?,Spinb?,Spinc?) * GammaId(Spind?,Spinb?) = Gamma(Mua,Spind,Spinc);
* Metric
    id Metric(Mua?,Mub?) * Gamma(Mua?,Spind?,Spinf?) = Gamma(Mub,Spind,Spinf);
    id Metric(Mua?,Mub?) * P(Mua?,Momd?) = P(Mub,Momd);
    id Metric(Mua?,Mub?) * Metric(Mub?,Mua?) = 4;
    id Metric(Mua?,Mub?) * Metric(Mub?,Muc?) = Metric(Mua,Muc);
    id Metric(Mua?,Mua?) = 4;
* standard Gamma algebra
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mua?,Spinc?,Spind?) = -GammaId(Mub,Mud);
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mud?,Spinc?,Spine?)*Gamma(Mua?,Spine?,Spinf?) = -2*Gamma(Mud,Spinb,Spinf);
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mud?,Spinc?,Spine?)*Gamma(Muf?,Spine?,Spinm?)*Gamma(Mua?,Spinm?,Spink?) = 4*Metric(Mud,Muf)*GammaId(Spinb,Spink);
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mud?,Spinc?,Spine?)*Gamma(Muf?,Spine?,Spinm?)*Gamma(Muk?,Spinm?,Spinl?)*Gamma(Mua?,Spinl?,Spinj?) = -2*Gamma(Muk,Spinb,Spinc)*Gamma(Muf,Spinc,Spinl)*Gamma(Mud,Spinl,Spinj);
* traces of Gamma
    id Gamma(Mua?,Spinb?,Spinb?) = 0;
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mub?,Spinc?,Spinb?) = 4*Metric(Mua,Mub);
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mub?,Spinc?,Spind?)*Gamma(Muc?,Spind?,Spinb?) = 0;
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mub?,Spinc?,Spind?)*Gamma(Muc?,Spind?,Spine?)*Gamma(Mud?,Spine?,Spinb?) 
        = 4*(Metric(Mua,Mub)*Metric(Muc,Mud) - Metric(Mua,Muc)*Metric(Mub,Mud)+ Metric(Mua,Mud)*Metric(Mub,Muc)) ;
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mub?,Spinc?,Spind?)*Gamma(Muc?,Spind?,Spine?)*Gamma(Mud?,Spine?,Spinf?)*Gamma(Mue?,Spinf?,Sping?)*Gamma(Muf?,Sping?,Spinb?)
        = 4*(Metric(Mua,Mub)*Metric(Muc,Mud)*Metric(Mue,Muf) 
           - Metric(Mua,Mub)*Metric(Muc,Mue)*Metric(Mud,Muf) 
           + Metric(Mua,Mub)*Metric(Muc,Muf)*Metric(Mud,Mue) 
           - Metric(Mua,Muc)*Metric(Mub,Mud)*Metric(Mue,Muf) 
           + Metric(Mua,Muc)*Metric(Mub,Mue)*Metric(Mud,Muf) 
           - Metric(Mua,Muc)*Metric(Mub,Muf)*Metric(Mud,Mue) 
           + Metric(Mua,Mud)*Metric(Mub,Muc)*Metric(Mue,Muf) 
           - Metric(Mua,Mud)*Metric(Mub,Mue)*Metric(Muc,Muf) 
           + Metric(Mua,Mud)*Metric(Mub,Muf)*Metric(Muc,Mue)
           - Metric(Mua,Mue)*Metric(Mub,Muc)*Metric(Mud,Muf)
           + Metric(Mua,Mue)*Metric(Mub,Mud)*Metric(Muc,Muf)
           - Metric(Mua,Mue)*Metric(Mub,Muf)*Metric(Muc,Mud)
           + Metric(Mua,Muf)*Metric(Mub,Muc)*Metric(Mud,Mue)
           - Metric(Mua,Muf)*Metric(Mub,Mud)*Metric(Muc,Mue)
           + Metric(Mua,Muf)*Metric(Mub,Mue)*Metric(Muc,Mud));
endrepeat;
"""

# TODO figure out why this does not work but manual dummies work
dirac_trick = """
repeat;
    id u(Spinc?,Momb?)*ubar(Spina?,Momb?) = Gamma(N5_?,Spinc,Spina) * P(N5_?,Momb) + GammaId(Spinc,Spina) * P(N5_?,Momb) * P(N5_?,Momb);
    id vbar(Spinc?,Momb?)*v(Spina?,Momb?) = Gamma(N7_?,Spinc,Spina) * P(N7_?,Momb) - GammaId(Spinc,Spina) * P(N7_?,Momb) * P(N7_?,Momb);
endrepeat;
"""

# TODO implement collecting of gammas and form calc solving of it


def get_gammas():
    return get_polarisation_sum() + get_dirac_trick() + gammas


def apply_gammas(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_gammas())


def get_polarisation_sum():
    polsum_feyn = """
    id epsstar(Muc?,Polb?,Moma?) * eps(Mul?,Pold?,Moma?) = -Metric(Muc,Mul);
    """
    polsum_feyn = """
    id epsstar(Muc?,Polb?,Moma?) * eps(Mul?,Pold?,Moma?) = -Metric(Mul,Muc);
    """
    return polsum_feyn


# TODO: look above
# def get_dirac_trick():
#    return dirac_trick


def get_dirac_trick(N=10):
    ret = ""
    for i in range(N):
        dummy = get_dummy_index()
        dirac_trick = f"""
    once u(Spinc?,Momb?)*ubar(Spina?,Momb?) = Gamma({dummy},Spinc,Spina) * P({dummy},Momb) + GammaId(Spinc,Spina) * P({dummy},Momb) * P({dummy},Momb);
    """
        ret += dirac_trick
    for i in range(N):
        dummy = get_dummy_index()
        dirac_trick = f"""
    once vbar(Spinc?,Momb?)*v(Spina?,Momb?) = Gamma({dummy},Spinc,Spina) * P({dummy},Momb) - GammaId(Spinc,Spina) * P({dummy},Momb) * P({dummy},Momb);
    """
        ret += dirac_trick
    return ret


def apply_dirac_trick(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_dirac_trick())


def apply_polarisation_sum(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_polarisation_sum())


# TODO: implement this use forms gamma algebra
def replace_indices_by_line():
    # return list of lines
    pass
