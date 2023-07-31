from feynpy.util import find_particle_in_model


def get_leg_math(fd,leg,model): # epsilons or u/v optionally also barred
    p = find_leg_in_model(fd,leg,model) 

    if p.spin == 3: 
        if leg.is_incoming():
            return f"Eps(Mu({p.particle.id}),Mom({p.particle.id}),Pol({p.particle.id}))"
        else:
            return f"Eps_star(Mu({p.particle.id}),Mom({p.particle.id}),Pol({p.particle.id}))"
    if p.spin == 2: 
        if not p.particle.is_anti():
            if leg.is_incoming():
                return f"U(Spin({p.particle.id}),Mom({p.particle.id}))"
            else:
                return f"U_bar(Spin({p.particle.id}),Mom({p.particle.id}))"
        else:
            if leg.is_incoming():
                return f"V(Spin({p.particle.id}),Mom({p.particle.id}))"
            else:
                return f"V_bar(Spin({p.particle.id}),Mom({p.particle.id}))"





def find_leg_in_model(fd,leg,model): # find leg in model
    assert leg in fd.legs
    return find_particle_in_model(leg,model)
