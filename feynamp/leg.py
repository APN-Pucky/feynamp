from feynamp.momentum import insert_momentum
from feynamp.util import find_particle_in_model


def get_leg_math_string(leg, fd, model):
    return get_leg_math(leg, fd, model)


def get_leg_math(fd, leg, model):  # epsilons or u/v optionally also barred
    p = find_leg_in_model(fd, leg, model)
    if p.particle.momentum is None or p.particle.momentum.name is None:
        raise ValueError("Momentum not set for particle")
    ret = ""
    mom = insert_momentum(p.particle.momentum.name)
    # give particles color vectors to sum over them in the end (or better average)
    # TODO this could be also done as incoming vs outcoming
    if p.color == 8:
        # if particle is a gluon give it a adjoint color function
        ret += f"VA(Glu{p.particle.id},{mom})*"
    if p.color == 3 or p.color == -3:
        # if particle is a quark give it a fundamental color function
        ret += f"VC(Color{p.particle.id},{mom})*"

    if p.spin == 3:
        if leg.is_incoming():
            ret += f"eps(Mu{p.particle.id},Pol{p.particle.id},{mom})"
        else:
            ret += f"eps_star(Mu{p.particle.id},Pol{p.particle.id},{mom})"
    if p.spin == 2:
        if not p.particle.is_anti():
            if leg.is_incoming():
                ret += f"u(Spin{p.particle.id},{mom})"
            else:
                ret += f"u_bar(Spin{p.particle.id},{mom})"
        else:
            if leg.is_incoming():
                ret += f"v(Spin{p.particle.id},{mom})"
            else:
                ret += f"v_bar(Spin{p.particle.id},{mom})"
    return ret


def find_leg_in_model(fd, leg, model):  # find leg in model
    assert leg in fd.legs
    return find_particle_in_model(leg, model)
