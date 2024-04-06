import re
from typing import List

import numpy as np
from feynml.connector import Connector
from feynml.id import generate_new_id
from feynml.leg import Leg
from feynml.propagator import Propagator
from feynml.vertex import Vertex

from feynamp.log import debug, info
from feynamp.momentum import insert_momentum
from feynamp.util import safe_index_replace


def insert_color_types(s: str):
    """ """
    # We use the non greedy .*? to match multiple occurances individually
    # s = re.sub(r"T\((.*?),(.*?),(.*?)\)", r"T(Glu\1,Color\2,Color\3)", s)
    # TODO do we want this? Or above
    s = re.sub(r"T\((.*?),(.*?),(.*?)\)", r"T(Color\2,Color\3,Glu\1)", s)
    s = re.sub(r"f\((.*?),(.*?),(.*?)\)", r"f(Glu\1,Glu\2,Glu\3)", s)
    s = re.sub(r"Identity\((.*?),(.*?)\)", r"Identity(Color\1,Color\2)", s)
    return s


def insert_lorentz_types(s: str, connections: List[Connector], vertex: Vertex):
    # We use the non greedy .*? to match multiple occurances individually
    s = re.sub(r"Gamma\((.*?),(.*?),(.*?)\)", r"Gamma(Mu\1,Spin\2,Spin\3)", s)
    s = re.sub(r"ProjP\((.*?),(.*?)\)", r"ProjP(Spin\1,Spin\2)", s)
    s = re.sub(r"ProjM\((.*?),(.*?)\)", r"ProjM(Spin\1,Spin\2)", s)
    s = re.sub(r"Metric\((.*?),(.*?)\)", r"Metric(Mu\1,Mu\2)", s)
    # use insert_momentum to replace the second argument to P
    matches = re.findall(r"P\((.*?),(.*?)\)", s)
    for g in matches:
        # find connection with g[1] id in connections
        repl = None
        for c in connections:
            if c.id == g[1] or "In" + str(c.id) == g[1] or "Out" + str(c.id) == g[1]:
                # TODO: is this correct?
                if c.goes_into(vertex):
                    repl = (
                        "(P(Mu" + g[0] + "," + insert_momentum(c.momentum.name) + "))"
                    )
                    break
                else:
                    repl = (
                        "(-P(Mu" + g[0] + "," + insert_momentum(c.momentum.name) + "))"
                    )
                    break
        if repl is None:
            raise Exception(
                f"Connection with id {g[1]} not found in connections {connections}"
            )
        s = s.replace("P(" + g[0] + "," + g[1] + ")", repl)
    return s


# def insert_index_types(s):
#    s = insert_color_types(s)
#    s = insert_lorentz_types(s)
#    return s


def get_vertex_math_string(fd, vertex, model):
    vv = get_vertex_math(fd, vertex, model)
    s = ""
    for v in vv:
        s += f"({v[0]})*({v[1]})*({v[2]}) + "
    return s[:-3]


def get_index_prefix(connection: Connector, vertex):
    if connection.goes_into(vertex):
        return "In"
    else:
        return "Out"


def get_vertex_math(fd, vertex, model, typed=True):  # TODO subst negative indices
    if not typed:
        raise NotImplementedError("Only typed vertices are supported")
    vv = fd.get_connections(vertex)
    v = find_vertex_in_model(fd, vertex, model)
    if v is None:
        raise Exception(f"Vertex {vertex} not found in model")
        # return None
    assert len(v.color) == len(v.lorentz)
    cret = []
    lret = []
    info(f"{v.color=}")
    info(f"{v.lorentz=}")
    for j in range(len(v.color)):
        col = v.color[j]
        nid = generate_new_id()
        col = safe_index_replace(col, str(-1), str(nid))
        for i, vv in enumerate(v.particles):
            if isinstance(v.connections[i], Leg):
                col = safe_index_replace(col, str(i + 1), str(v.connections[i].id))
            elif isinstance(v.connections[i], Propagator):
                col = safe_index_replace(
                    col,
                    str(i + 1),
                    get_index_prefix(v.connections[i], vertex)
                    + str(v.connections[i].id),
                )
            else:
                raise Exception(
                    f"Connection {v.connections[i]} not a leg or propagator"
                )
        if typed:
            col = insert_color_types(col)
        cret.append(col)
    for k in range(len(v.lorentz)):
        lor = v.lorentz[k].structure
        debug(f"{lor=}")
        nid = generate_new_id()
        lor = safe_index_replace(lor, str(-1), str(nid))
        for i, vv in enumerate(v.particles):
            if isinstance(v.connections[i], Leg):
                lor = safe_index_replace(lor, str(i + 1), str(v.connections[i].id))
            elif isinstance(v.connections[i], Propagator):
                lor = safe_index_replace(
                    lor,
                    str(i + 1),
                    get_index_prefix(v.connections[i], vertex)
                    + str(v.connections[i].id),
                )
            else:
                raise Exception(
                    f"Connection {v.connections[i]} not a leg or propagator"
                )
        if typed:
            lor = insert_lorentz_types(lor, v.connections, vertex)
        lret.append(lor)
    vertex_math = []
    for k, v in v.couplings.items():
        vertex_math.append((v.value, cret[k[0]], lret[k[1]]))
    debug(f"{vertex_math=}")
    return vertex_math


def find_vertex_in_model(fd, vertex, model):
    """
    Finds the model vertex corresponding to the given FeynmanDiagram vertex

    Note: Sorting is to check for the correct particles in a vertex given they can be in any order and have duplicates
    """
    assert vertex in fd.vertices
    cons = np.array(fd.get_connections(vertex))
    # debug(f"{cons=}")
    pdg_ids_list = []

    # correct for incoming vs outgoing fermion struct
    for c in cons:
        p = c.pdgid
        if c.is_any_fermion():
            if c.goes_into(vertex):
                p = -p
        pdg_ids_list += [p]
    pdg_ids_array = np.array(pdg_ids_list)

    sort_mask = np.argsort(pdg_ids_array)
    particles = pdg_ids_array[sort_mask]
    scons = cons[sort_mask]
    # debug(f"{scons=}")
    ret = None
    for v in model.vertices:
        if len(v.particles) != len(particles):
            continue
        model_particle_ids = np.array([p.pdg_code for p in v.particles])
        model_sort_mask = np.argsort(model_particle_ids)
        # By sorting based on the indices we reproduce the order of the particles in the vertex
        inverted_model_sort_mask = np.argsort(model_sort_mask)
        sorted_model_particle_ids = model_particle_ids[model_sort_mask]
        if np.array_equal(sorted_model_particle_ids, particles):
            vc = []
            for i, _ in enumerate(model_particle_ids):
                con = scons[inverted_model_sort_mask[i]]
                vc.append(con)
            v.connections = vc
            ret = v
            break

    # Make sure all connections are in the vertex
    for c in cons:
        assert c in ret.connections
    # debug(f"{ret=}")
    if ret is None:
        raise Exception(
            f"Vertex {vertex} with cons {cons} not found in model\n{pdg_ids_list=}"
        )
    return ret
