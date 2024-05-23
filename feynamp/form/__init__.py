import os
import re
from typing import List

import form
import sympy
from feynml.feynmandiagram import FeynmanDiagram
from feynml.leg import Leg
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
from feynamp.leg import get_leg_momentum
from feynamp.log import debug

# TODO compute squared  functino which coutns legs!!"!!!" and picks right mandelstamm,s


def compute_squared(fds: List[FeynmanDiagram], fm: FeynModel, tag=False):
    return compute_squared_correlated(
        fds, fm, colore_leg1=None, color_leg2=None, tag=tag
    )


def compute_squared_correlated(
    fds: List[FeynmanDiagram],
    fm: FeynModel,
    color_leg1: Leg,
    color_leg2: Leg,
    tag=False,
):
    assert len(fds) > 0, "No FeynmanDiagrams to compute"
    assert not (
        color_leg1 is not None and color_leg2 is None
    ), "If you want to color correlate legs, you need to provide both"
    assert not (
        color_leg1 is None and color_leg2 is not None
    ), "If you want to color correlate legs, you need to provide both"
    dims = fds[0].get_externals_size()
    for fd in fds:
        assert (
            dims == fd.get_externals_size()
        ), "All FeynmanDiagrams must have the same external legs"
    s2 = amplitude.square_parallel(fds, fm, tag=tag)
    debug(f"{s2=}")

    s2 = apply_color_parallel(
        s2, get_leg_momentum(color_leg1), get_leg_momentum(color_leg2)
    )

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
    fs += get_gammas(fds, fm)
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
