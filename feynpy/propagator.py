from feynml.id import generate_new_id

from feynpy.momentum import insert_momentum
from feynpy.util import find_particle_in_model


def get_propagator_math_string(fd, prop, model):
    return get_propagator_math(fd, prop, model)


def get_propagator_math(fd, prop, model):
    # find the particle in the model
    p = find_propagator_in_model(fd, prop, model)
    if p.particle.momentum is None or p.particle.momentum.name is None:
        raise ValueError("Momentum not set for particle")
    mom = insert_momentum(p.particle.momentum.name)
    # if boson just 1/(p^2-m^2)
    if p.spin == 3:
        # nid = generate_new_id()
        # TODO treate denominators differently for loops etc?
        return f"Denom({mom},{p.mass.name})"
    if p.spin == 2:  # TODO handle plus minus mass for fermions
        nid = generate_new_id()
        return f"(P(Mu{nid},{mom})*Gamma(Mu{nid},Spin{p.particle.source},Spin{p.particle.target}) + {p.mass.name}*GammaId(Spin{p.particle.source},Spin{p.particle.target}))*Denom({mom},{p.mass.name})"
    raise ValueError("Spin not set for particle")


def find_propagator_in_model(fd, prop, model):
    assert prop in fd.propagators
    return find_particle_in_model(prop, model)
