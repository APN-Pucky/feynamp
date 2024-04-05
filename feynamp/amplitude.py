import re
from typing import List

from feynml.feynmandiagram import FeynmanDiagram
from feynmodel.feyn_model import FeynModel

from feynamp.leg import get_leg_math_string
from feynamp.log import debug
from feynamp.propagator import get_propagator_math_string
from feynamp.vertex import get_vertex_math_string


def complex_conjugate(s: str):
    """"""
    # return s
    return re.sub(r"complex\((.*?),(.*?)\)", r"complex(\1,-\2)", s)


def feynman_diagram_to_string(feynman_diagram, feyn_model):
    fd = feynman_diagram
    # TODO use these informations from qgraf, might need to be implemented in feynml
    # sign = fd.get_sign()
    # symmetry = fd.get_symmetry_factor()
    vm = []
    lm = []
    pm = []
    for v in fd.vertices:
        vm.append("(" + get_vertex_math_string(fd, v, feyn_model) + ")")
    for l in fd.legs:
        lm.append("(" + get_leg_math_string(fd, l, feyn_model) + ")")
    for p in fd.propagators:
        pm.append("(" + get_propagator_math_string(fd, p, feyn_model) + ")")

    ret = ""
    if len(vm) > 0:
        ret += f"{' * '.join(vm)} * "
    if len(lm) > 0:
        ret += f"{' * '.join(lm)} * "
    if len(pm) > 0:
        ret += f"{' * '.join(pm)} * "
    debug(f"{ret=}")
    return ret[0:-3]


def multiply(
    lst_fd1: List[FeynmanDiagram], lst_fd2: List[FeynmanDiagram], feyn_model: FeynModel
):
    # TODO should this care about fermion lines!?
    s = ""
    lst_fd1 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd1]
    lst_fd2 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd2]
    for fd1 in lst_fd1:
        for fd2 in lst_fd2:
            s += f"({fd1})*({fd2}) + "
    return s[:-3]


def square(lst_fd: List[FeynmanDiagram], feyn_model: FeynModel, tag=False) -> str:
    """
    Squares the list of feynman diagrams taking the fermion sign into account.
    """
    return " + ".join(square_parallel(lst_fd, feyn_model, tag))


def square_parallel(
    lst_fd: List[FeynmanDiagram], feyn_model: FeynModel, tag=False
) -> List[str]:
    """
    Squares the list of feynman diagrams taking the fermion sign into account.
    """
    dims = lst_fd[0].get_externals_size()
    for fd in lst_fd:
        assert (
            dims == fd.get_externals_size()
        ), "All FeynmanDiagrams must have the same external legs"
    # TODO handle relative fermion sign (also majorana!) https://cds.cern.ch/record/238903/files/th-6549-92.pdf
    # TODO also multiply by the symmetry factor from qgraf
    # return multiply(lst_fd,[l.conjugated() for l in lst_fd],feyn_model)
    s = ""
    lst_fd1 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd]
    lst_fd2 = [
        complex_conjugate(feynman_diagram_to_string(l.conjugated(), feyn_model))
        for l in lst_fd
    ]
    debug(f"{lst_fd1=}")
    debug(f"{lst_fd2=}")
    ret_lst = []
    # TODO this could also be done in multiply by comparing the diagrams
    for i in range(len(lst_fd1)):
        # TODO reenable loop from i
        for j in range(len(lst_fd2)):
            sfd1 = lst_fd1[i]
            sfd2 = lst_fd2[j]
            # TODO reenable
            # if i == j:
            ttag = ""
            if tag:
                ttag = f"*fd{lst_fd[i].id}*fd{lst_fd[j].id}*fd{lst_fd[i].id}fd{lst_fd[j].id}"
            ferm_fac = lst_fd[i].get_fermion_factor(lst_fd[j])
            debug(f"{ferm_fac=}")
            ret_lst.append(f"({sfd1})*({sfd2}){ttag}*{ferm_fac}")
            # TODO reenable
            # elif i < j:
            #    ttag = ""
            #    if tag:
            #        ttag = f"*fd{lst_fd[i].id}*fd{lst_fd[j].id}*fd{lst_fd[i].id}fd{lst_fd[j].id}"
            #    ferm_fac = lst_fd[i].get_fermion_factor(lst_fd[j])
            #    debug(f"{ferm_fac=}")
            #    ret_lst.append(
            #        f"2*(+{sfd1})*({sfd2}){ttag}*{ferm_fac}"
            #    )  # TODO this needs Re!
    return ret_lst


def add(lst_fd, feyn_model):
    lst_fd1 = [feynman_diagram_to_string(l, feyn_model) for l in lst_fd]
    ferm_facs = [lst_fd[0].get_fermion_factor(l) for l in lst_fd]
    s = ""
    for i in range(len(lst_fd1)):
        s += f"({lst_fd1[i]})*{ferm_facs[i]} + "
    return s[:-3]
