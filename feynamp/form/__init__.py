import os
import re
from typing import List

import form
import sympy
from feynml.feynmandiagram import FeynmanDiagram
from feynmodel.feyn_model import FeynModel

from feynamp import get_color_average, get_spin_average
import feynamp.amplitude as amplitude
from feynamp.log import debug
from feynamp.form.color import get_color
from feynamp.form.lorentz import get_gammas, get_polarisation_sums
from feynamp.form.momentum import get_kinematics, get_onshell, get_mandelstamm, apply, apply_den

# TODO compute squared  functino which coutns legs!!"!!!" and picks right mandelstamm,s


def compute_squared(fds: List[FeynmanDiagram], fm: FeynModel):
    dims = fds[0].get_externals_size()
    for fd in fds:
        assert (
            dims == fd.get_externals_size()
        ), "All FeynmanDiagrams must have the same external legs"
    s2 = amplitude.square(fds, fm, tag=False)
    debug(f"{s2=}")
    fs = ""
    fs += get_polarisation_sums(fds, fm)
    fs += get_gammas()
    fs += get_color()
    fs += get_kinematics()
    fs += get_onshell(fds, fm)
    fs += get_mandelstamm(fds, fm)

    rs = apply(s2, fs)
    debug(f"{rs=}")

    rr = apply_den(
        rs,
        get_onshell(fds, fm)
        + get_mandelstamm(fds, fm),
    )
    debug(f"{rr=}")

    ret = sympy.simplify(
        sympy.parse_expr(
            rr.replace("Mom_", "")
            .replace(".", "_")
            .replace("^", "**")
            .replace("mss", "s")
            .replace("msu", "u")
            .replace("mst", "t")
        )
        * sympy.parse_expr(
            "*".join([*get_color_average(fds), *get_spin_average(fds)])
        )
    )
    return ret
