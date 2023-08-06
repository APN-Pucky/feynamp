from sympy.parsing.sympy_parser import parse_expr

from feynpy.leg import get_leg_math_string
from feynpy.lorentz import gamma
from feynpy.propagator import get_propagator_math_string
from feynpy.vertex import get_vertex_math_string


def feynman_diagram_to_string(feynman_diagram, feyn_model):
    fd = feynman_diagram
    vm = []
    lm = []
    pm = []
    for v in fd.vertices:
        vm.append(get_vertex_math_string(fd, v, feyn_model))
    for l in fd.legs:
        lm.append(get_leg_math_string(fd, l, feyn_model))
    for p in fd.propagators:
        pm.append(get_propagator_math_string(fd, p, feyn_model))
    return f"{' * '.join(vm)} * {' * '.join(lm)} * {' * '.join(pm)}"


def string_to_sympy(s):
    s = s.replace(
        "Gamma", "gamma"
    )  # we keep string like the ufo, but want lowercase gamma for sympy
    s = s.replace("Metric", "metric")
    s = s.replace("complex(0,1)", "I")  # sympy uses I for imaginary unit
    return parse_expr(s, local_dict={"gamma": gamma})
    # return parse_expr(s, evaluate=False)


def multiply(lst_fd1, lst_fd2, feyn_model):
    s = ""
    lst_fd1 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd1]
    lst_fd2 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd2]
    for fd1 in lst_fd1:
        for fd2 in lst_fd2:
            s += f"({fd1})*({fd2}) + "
    return s[:-3]


def square(lst_fd, feyn_model, tag=True):
    # TODO handle relative fermion sign (also majorana!) https://cds.cern.ch/record/238903/files/th-6549-92.pdf
    # return multiply(lst_fd,[l.conjugated() for l in lst_fd],feyn_model)
    s = ""
    lst_fd1 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd]
    lst_fd2 = [feynman_diagram_to_string(l.conjugated(), feyn_model) for l in lst_fd]
    # TODO this could also be done in multiply by comparing the diagrams
    for i in range(len(lst_fd1)):
        for j in range(i, len(lst_fd2)):
            sfd1 = lst_fd1[i]
            sfd2 = lst_fd2[j]
            if i == j:
                ttag = ""
                if tag:
                    ttag = f"*fd{lst_fd[i].id}*fd{lst_fd[i].id}fd{lst_fd[i].id}"
                s += f"({sfd1})*({sfd2}){ttag} + "
            elif i < j:
                ttag = ""
                if tag:
                    ttag = f"*fd{lst_fd[i].id}*fd{lst_fd[j].id}*fd{lst_fd[i].id}fd{lst_fd[j].id}"
                s += f"2*(+{sfd1})*({sfd2}){ttag} + "  # TODO this needs Re
    return s[:-3]
