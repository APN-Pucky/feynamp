from feynpy.util import find_particle_in_model


def get_leg_math_string(leg, fd, model):
    return get_leg_math(leg, fd, model)


def get_leg_math(fd, leg, model, typed=True):  # epsilons or u/v optionally also barred
    p = find_leg_in_model(fd, leg, model)

    if p.spin == 3:
        if leg.is_incoming():
            if typed:
                return f"Eps(Mu({p.particle.id}),Mom({p.particle.id}),Pol({p.particle.id}))"
            else:
                return f"Eps_star({p.particle.id},{p.particle.id},{p.particle.id})"
        else:
            if typed:
                return f"Eps_star(Mu({p.particle.id}),Mom({p.particle.id}),Pol({p.particle.id}))"
            else:
                return f"Eps_star({p.particle.id},{p.particle.id},{p.particle.id})"
    if p.spin == 2:
        if not p.particle.is_anti():
            if leg.is_incoming():
                if typed:
                    return f"U(Spin({p.particle.id}),Mom({p.particle.id}))"
                else:
                    return f"U({p.particle.id},{p.particle.id})"
            else:
                if typed:
                    return f"U_bar(Spin({p.particle.id}),Mom({p.particle.id}))"
                else:
                    return f"U_bar({p.particle.id},{p.particle.id})"
        else:
            if leg.is_incoming():
                if typed:
                    return f"V(Spin({p.particle.id}),Mom({p.particle.id}))"
                else:
                    return f"V({p.particle.id},{p.particle.id})"
            else:
                if typed:
                    return f"V_bar(Spin({p.particle.id}),Mom({p.particle.id}))"
                else:
                    return f"V_bar({p.particle.id},{p.particle.id})"


def find_leg_in_model(fd, leg, model):  # find leg in model
    assert leg in fd.legs
    return find_particle_in_model(leg, model)
