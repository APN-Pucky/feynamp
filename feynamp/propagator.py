from feynml.id import generate_new_id

from feynamp.momentum import insert_mass, insert_momentum
from feynamp.util import find_particle_in_model


def get_propagator_math_string(fd, prop, model):
    return get_propagator_math(fd, prop, model)


def get_propagator_math(fd, prop, model):
    # find the particle in the model
    p = find_propagator_in_model(fd, prop, model)
    if p.particle.momentum is None or p.particle.momentum.name is None:
        raise ValueError("Momentum not set for particle")
    mom = insert_momentum(p.particle.momentum.name)
    mass = insert_mass(p.mass.name)
    ret = ""
    if p.color == 3 or p.color == -3:
        ret += f"df(ColorIn{p.particle.id},ColorOut{p.particle.id})*"
    if p.color == 8:
        ret += f"da(GluIn{p.particle.id},GluOut{p.particle.id})*"
    # if boson just 1/(p^2-m^2)
    if p.spin == 3:
        # nid = generate_new_id()
        # TODO treate denominators differently for loops etc?
        ret += f"(-1)*complex(0,1)*Metric(MuIn{p.particle.id},MuOut{p.particle.id})*"
    elif p.spin == 2:  # TODO handle plus minus mass for fermions
        nid = generate_new_id()
        ret += f"complex(0,1)*(P(Mu{nid},{mom})*Gamma(Mu{nid},SpinIn{p.particle.id},SpinOut{p.particle.id})"
        ret += f" + {mass}*GammaId(SpinIn{p.particle.id},SpinOut{p.particle.id}))*"
    else:
        raise ValueError("Spin not set for particle")
    return ret + f"Denom({mom},{mass})"


def find_propagator_in_model(fd, prop, model):
    assert prop in fd.propagators
    return find_particle_in_model(prop, model)
