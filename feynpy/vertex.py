import re

import numpy as np
from feynml.id import generate_new_id
from feynml.leg import Leg
from feynml.propagator import Propagator

from feynpy.momentum import insert_momentum
from feynpy.util import safe_index_replace


def insert_color_types(s):
    s = re.sub(r"T\((.*),(.*),(.*)\)", r"T(Glu\1,Col\2,Col\3)", s)
    s = re.sub(r"f\((.*),(.*),(.*)\)", r"f(Glu\1,Glu\2,Glu\3)", s)
    return s


def insert_lorentz_types(s):
    s = re.sub(r"Gamma\((.*),(.*),(.*)\)", r"Gamma(Mu\1,Spin\2,Spin\3)", s)
    s = re.sub(r"ProjP\((.*),(.*)\)", r"ProjP(Spin\1,Spin\2)", s)
    s = re.sub(r"ProjM\((.*),(.*)\)", r"ProjM(Spin\1,Spin\2)", s)
    s = re.sub(r"Metric\((.*),(.*)\)", r"Metric(Mu\1,Mu\2)", s)
    # use insert_momentum to replace the second argument to P
    m = re.match(r"P\((.*),(.*)\)", s)
    if m:
        for g in m.groups():
            s.replace(
                "P(" + g[0] + "," + g[1] + ")",
                "P(Mu" + g[0] + "," + insert_momentum(g[1]) + ")",
            )
    return s


def insert_index_types(s):
    s = insert_color_types(s)
    s = insert_lorentz_types(s)
    return s


def get_vertex_math_string(fd, vertex, model):
    vv = get_vertex_math(fd, vertex, model)
    s = ""
    for v in vv:
        s += f"{v[0]}*{v[1]}*{v[2]} + "
    return s[:-3]


def get_vertex_math(fd, vertex, model, typed=True):  # TODO subst negative indices
    vv = fd.get_connections(vertex)
    v = find_vertex_in_model(fd, vertex, model)
    if v is None:
        return None
    assert len(v.color) == len(v.lorentz)
    cret = []
    lret = []
    for j in range(len(v.color)):
        col = v.color[j]
        nid = generate_new_id()
        col = safe_index_replace(col, str(-1), str(nid))
        for i, vv in enumerate(v.particles):
            col = safe_index_replace(col, str(i + 1), str(v.connections[i].id))
        if typed:
            col = insert_color_types(col)
        cret.append(col)
    for k in range(len(v.lorentz)):
        lor = v.lorentz[j].structure
        nid = generate_new_id()
        lor = safe_index_replace(lor, str(-1), str(nid))
        for i, vv in enumerate(v.particles):
            lor = safe_index_replace(lor, str(i + 1), str(v.connections[i].id))
        if typed:
            lor = insert_lorentz_types(lor)
        lret.append(lor)
    ret = []
    for k, v in v.couplings.items():
        ret.append((v.value, cret[k[0]], lret[k[1]]))
    return ret


def find_vertex_in_model(fd, vertex, model):
    # TODO handle multiple vertices
    assert vertex in fd.vertices
    cons = np.array(fd.get_connections(vertex))
    aa = []
    for c in cons:
        p = c.pdgid
        # correct for incoming vs outgoing fermion struct
        if c.is_any_fermion():
            if c.goes_into(vertex):
                p = -p

        aa += [p]
    cpd = np.array(aa)

    cmask = np.argsort(cpd)
    particles = cpd[cmask]
    scons = cons[cmask]
    for v in model.vertices:
        if len(v.particles) != len(particles):
            continue
        pp = np.array([p.pdg_code for p in v.particles])
        smp = sorted(pp)
        if np.array_equal(smp, particles):
            vc = []
            for i, ps in enumerate(pp):
                vc.append(scons[smp.index(ps)])
            v.connections = vc
            return v
    return None
