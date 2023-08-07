from feynamp.leg import get_leg_math_string
from feynamp.propagator import get_propagator_math_string
from feynamp.vertex import get_vertex_math_string


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


def multiply(lst_fd1, lst_fd2, feyn_model):
    # TODO should this care about fermion lines!?
    s = ""
    lst_fd1 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd1]
    lst_fd2 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd2]
    for fd1 in lst_fd1:
        for fd2 in lst_fd2:
            s += f"({fd1})*({fd2}) + "
    return s[:-3]


def square(lst_fd, feyn_model, tag=False):
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
                ferm_fac = lst_fd[i].get_fermion_factor(lst_fd[j])
                s += f"2*(+{sfd1})*({sfd2}){ttag}*{ferm_fac} + "  # TODO this needs Re!
    return s[:-3]


def add(lst_fd, feyn_model):
    lst_fd1 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd]
    ferm_facs = [lst_fd[0].get_fermion_factor(l) for l in lst_fd]
    s = ""
    for i in range(len(lst_fd1)):
        s += f"({lst_fd1[i]})*{ferm_facs[i]} + "
    return s[:-3]
