from typing import List

from feynml.feynmandiagram import FeynmanDiagram
from feynml.leg import Leg
from feynmodel.feyn_model import FeynModel

from feynamp.form.form import get_dummy_index, init, run, string_to_form
from feynamp.leg import find_leg_in_model
from feynamp.log import debug
from feynamp.momentum import insert_momentum
from feynamp.util import is_mass_zero

gammas = """
repeat;
* identity
    id Gamma(Mua?,Spinb?,Spinc?) * GammaId(Spinc?,Spind?) = Gamma(Mua,Spinb,Spind);
    id Gamma(Mua?,Spinb?,Spinc?) * GammaId(Spind?,Spinb?) = Gamma(Mua,Spind,Spinc);
    id GammaId(Spina?,Spinb?) * GammaId(Spinb?,Spind?) = GammaId(Spina,Spind);
    id GammaId(Spina?,Spina?) = 4;
* Metric
    id Metric(Mua?,Mub?) * Gamma(Mua?,Spind?,Spinf?) = Gamma(Mub,Spind,Spinf);
    id Metric(Mua?,Mub?) * P(Mua?,Momd?) = P(Mub,Momd);
    id Metric(Mua?,Mub?) * Metric(Mub?,Mua?) = 4;
    id Metric(Mua?,Mub?) * Metric(Mub?,Muc?) = Metric(Mua,Muc);
    id Metric(Mua?,Mua?) = 4;
* standard Gamma algebra
    id Gamma(Mua?,Spinb?,Spinc?)*Gamma(Mua?,Spinc?,Spind?) = -GammaId(Spinb,Spind);
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


# TODO implement collecting of gammas and form calc solving of it
# idea: GammaCollect(1,spin1,spin2, mu,nu,...) run through and then apply to expression
gamma_collect = """
#do i = 1, 10
once Gamma(Mux?,Spin1?,Spin2?) = GammaCollect(`i',Spin1,Spin2,Mux);
repeat;
id GammaCollect(`i',Spin1?,Spin2?,?mus)*Gamma(Mux?, Spin2?,Spin3?) = GammaCollect(`i',Spin1,Spin3, ?mus, Mux);
id GammaCollect(`i',Spin1?,Spin2?,?mus)*GammaId(Spin2?,Spin3?) = GammaCollect(`i',Spin1,Spin3, ?mus)*gi_(`i');
endrepeat;
id GammaCollect(`i',Spin1?,Spin1?,?mus) = g_(`i',?mus);
trace4, `i';
#enddo
#do i = 11, 21
once GammaId(Spin1?,Spin2?) = GammaIdCollect(`i',Spin1,Spin2);
repeat;
id GammaIdCollect(`i',Spin1?,Spin2?)*GammaId(Spin2?,Spin3?) = GammaIdCollect(`i',Spin1,Spin3)*gi_(`i');
endrepeat;
id GammaIdCollect(`i',Spin1?,Spin1?) = 1;
trace4, `i';
#enddo
"""

metrics = """
* Simplify metrics, here let form handle it 
repeat;
id Metric(Mua?,Mub?) = d_(Mua,Mub);
endrepeat;
"""


def get_metrics():
    return metrics


def apply_metrics(string_expr):
    s = string_to_form(string_expr)
    from .color import color_init

    return run(init + color_init + f"Local TMP = {s};" + get_metrics())


def get_gammas():
    return new_get_gammas()


def new_get_gammas():
    return get_dirac_trick() + get_metrics() + gamma_collect


def get_gammas_v1():
    return get_polarisation_sum_v1() + get_dirac_trick() + gammas


def apply_gammas(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_gammas())


def get_orthogonal_polarisation_momentum(
    leg: Leg, fds: List[FeynmanDiagram], model: FeynModel
):
    for fd in [fds[0]]:
        for fleg in fd.legs:
            p = find_leg_in_model(fd, fleg, model)
            if (
                is_mass_zero(p)
                and fleg != leg
                and leg.is_incoming() == fleg.is_incoming()
            ):
                mom = insert_momentum(fleg.momentum.name)
                return mom
    raise ValueError("No orthogonal momentum found")


def get_polarisation_sums(fds: List[FeynmanDiagram], model: FeynModel):
    pol_sums = ""
    # TODO might want to loop over all fds?
    for fd in [fds[0]]:
        for l in fd.legs:
            p = find_leg_in_model(fd, l, model)
            mom = insert_momentum(l.momentum.name)
            # mass = insert_mass(string_to_form(p.mass.name))
            if p.spin == 3:
                if is_mass_zero(p):
                    if p.color == 8:
                        mom_n = get_orthogonal_polarisation_momentum(l, fds, model)
                        pol_sums += get_polarisation_sum_physical(mom, mom_n)
                        # pol_sums += get_polarisation_sum_feynman(mom)
                        debug(f"mom: {mom}, mom_n: {mom_n}")
                    elif p.color == 1:
                        pol_sums += get_polarisation_sum_feynman(mom)
                else:
                    pol_sums += get_polarisation_sum_massive(mom)
    debug(f"pol_sums: {pol_sums}")
    return pol_sums


def get_polarisation_sum_massive(mom_a):
    pol_sum = f"""
    id epsstar(Muc?,Polb?,{mom_a}) * eps(Mul?,Pold?,{mom_a}) = -Metric(Mul,Muc) + (P(Mul,{mom_a})*P(Muc,{mom_a}))*Den({mom_a}.{mom_a});
    """
    return pol_sum


def get_polarisation_sum_feynman(mom_a):
    pol_sum = f"""
    id epsstar(Muc?,Polb?,{mom_a}) * eps(Mul?,Pold?,{mom_a}) = -Metric(Mul,Muc);
    """
    return pol_sum


def get_polarisation_sum_physical(mom_a, mom_b):
    pol_sum = f"""
    id epsstar(Muc?,Polb?,{mom_a}) * eps(Mul?,Pold?,{mom_a}) = -Metric(Muc,Mul) 
    + (P(Muc,{mom_a})*P(Mul,{mom_b}) +  P(Mul,{mom_a})*P(Muc,{mom_b}))*Den({mom_b}.{mom_a}) 
    - P(Muc,{mom_a})*P(Mul,{mom_a})*({mom_b}.{mom_b})*Den({mom_b}.{mom_a})*Den({mom_b}.{mom_a});
    """
    return pol_sum


def get_polarisation_sum(mom_a, mom_b=None):
    if mom_b is None:
        # massive case
        return get_polarisation_sum_massive(mom_a)
    elif mom_b == 0:
        # photon massless case
        return get_polarisation_sum_feynman(mom_a)
    else:
        # gluon massless case
        return get_polarisation_sum_physical(mom_a, mom_b)


def get_polarisation_sum_v1():
    polsum_feyn = """
    id epsstar(Muc?,Polb?,Moma?) * eps(Mul?,Pold?,Moma?) = -Metric(Muc,Mul);
    """
    polsum_feyn = """
    id epsstar(Muc?,Polb?,Moma?) * eps(Mul?,Pold?,Moma?) = -Metric(Mul,Muc);
    """
    return polsum_feyn


# TODO figure out why this does not work but manual dummies work
dirac_trick = """
repeat;
    id u(Spinc?,Momb?)*ubar(Spina?,Momb?) = Gamma(N5_?,Spinc,Spina) * P(N5_?,Momb) + GammaId(Spinc,Spina) * P(N5_?,Momb) * P(N5_?,Momb);
    id vbar(Spinc?,Momb?)*v(Spina?,Momb?) = Gamma(N7_?,Spinc,Spina) * P(N7_?,Momb) - GammaId(Spinc,Spina) * P(N7_?,Momb) * P(N7_?,Momb);
endrepeat;
"""


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


# TODO: look above
# def get_dirac_trick():
#    return dirac_trick


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
