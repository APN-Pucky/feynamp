# Utilities for feynamp
from feynmodel.particle import Particle


def safe_index_replace(string, old, new):
    string = string.replace("," + old + ",", "," + new + ",")
    string = string.replace("(" + old + ",", "(" + new + ",")
    string = string.replace("," + old + ")", "," + new + ")")
    return string


def find_particle_in_model(particle, model):
    for pp in model.particles:
        if pp.pdg_code == particle.pdgid:
            pp.particle = particle
            return pp
    return None


def is_mass_zero(p: Particle):
    return p.mass.name == "ZERO" or float(p.mass.value) == 0.0
