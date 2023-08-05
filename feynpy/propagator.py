from feynml.id import generate_new_id

from feynpy.util import find_particle_in_model


def get_propagator_math_string(fd, prop, model):
    return get_propagator_math(fd, prop, model)


def get_propagator_math(fd, prop, model):
    # find the particle in the model
    p = find_propagator_in_model(fd, prop, model)
    # if boson just 1/(p^2-m^2)
    if p.spin == 3:
        nid = generate_new_id()
        # TODO treate denominators differently for loops etc?
        return f"Denom({p.particle.id},{p.mass.name})"
    if p.spin == 2:  # TODO handle plus minus mass for fermions
        nid = generate_new_id()
        return f"(P({nid},{p,particle.id})*Gamma({nid},{p.particle.source},{p.particle.target}) + {p.mass.name})*Denom({p.particle.id},{p.mass.name})"


def find_propagator_in_model(fd, prop, model):
    assert prop in fd.propagators
    return find_particle_in_model(prop, model)
