import os
import re
from typing import List

import form
import sympy
from feynml.feynmandiagram import FeynmanDiagram
from feynmodel.feyn_model import FeynModel

import feynamp.amplitude as amplitude
from feynamp import get_color_average, get_spin_average
from feynamp.form.color import apply_color_parallel
from feynamp.form.lorentz import get_gammas, get_metrics, get_polarisation_sums
from feynamp.form.momentum import (
    apply,
    apply_den,
    apply_parallel,
    get_kinematics,
    get_mandelstamm,
    get_onshell,
)
from feynamp.log import debug

# TODO compute squared  functino which coutns legs!!"!!!" and picks right mandelstamm,s


def compute_squared(fds: List[FeynmanDiagram], fm: FeynModel, tag=False):
    dims = fds[0].get_externals_size()
    for fd in fds:
        assert (
            dims == fd.get_externals_size()
        ), "All FeynmanDiagrams must have the same external legs"
    s2 = amplitude.square_parallel(fds, fm, tag=tag)
    debug(f"{s2=}")

    s2 = apply_color_parallel(s2)

    fs = ""
    fs += get_metrics()
    # fs += get_color()
    # fs += get_kinematics()
    # fs += get_onshell(fds, fm)
    # fs += get_mandelstamm(fds, fm)
    # This is where it gets expensive
    fs += get_polarisation_sums(fds, fm)
    fs += get_kinematics()
    fs += get_onshell(fds, fm)
    fs += get_gammas()
    fs += get_kinematics()
    fs += get_onshell(fds, fm)
    fs += get_mandelstamm(fds, fm)

    rs = apply_parallel(s2, fs)
    rs = " + ".join([f"({r})" for r in rs])
    debug(f"{rs=}")

    rr = apply_den(
        rs,
        get_onshell(fds, fm) + get_mandelstamm(fds, fm),
    )
    debug(f"{rr=}")

    ret = sympy.parse_expr(
        rr.replace("Mom_", "")
        .replace(".", "_")
        .replace("^", "**")
        .replace("mss", "s")
        .replace("msu", "u")
        .replace("mst", "t")
    ) * sympy.parse_expr("*".join([*get_color_average(fds), *get_spin_average(fds)]))
    return ret
