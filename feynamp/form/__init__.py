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
from feynamp.form.form import apply_parallel_v3
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


# def compute_squared(fds: List[FeynmanDiagram], fm: FeynModel, tag=False):
#    return compute_squared_correlated(
#        fds, fm, colorcorrelated=False, tag=tag
#    )


def compute_squared(
    fds: List[FeynmanDiagram],
    fm: FeynModel,
    colorcorrelated=False,
    spincorrelated=False,
    tag=False,
    drop_ms_prefix=False,
    optimize=False,
    only_result=True,
    re_for_interference=True,
):
    if optimize:
        assert not only_result
    if spincorrelated:
        assert not optimize, "Does not work with Tensors..."  # TODO optimize this
    assert len(fds) > 0, "No FeynmanDiagrams to compute"
    dims = fds[0].get_externals_size()
    for fd in fds:
        assert (
            dims == fd.get_externals_size()
        ), "All FeynmanDiagrams must have the same external legs"
    s2 = amplitude.square_parallel(
        fds, fm, tag=tag, prefactor=True, re_for_interference=re_for_interference
    )
    # debug(f"{s2=}")

    s2 = apply_color_parallel(
        s2, fds=fds, legs=fds[0].legs, model=fm, colorcorrelation=colorcorrelated
    )

    fs = ""
    fs += get_metrics()
    # fs += get_color()
    # fs += get_kinematics()
    # fs += get_onshell(fds, fm)
    # fs += get_mandelstamm(fds, fm)
    # This is where it gets expensive
    fs += get_polarisation_sums(fds, fm, spincorrelated=spincorrelated)
    fs += get_kinematics()
    fs += get_onshell(fds, fm)
    fs += get_gammas(fds, fm)
    fs += get_kinematics()
    fs += get_onshell(fds, fm)
    fs += get_mandelstamm(fds, fm)

    rs = apply_parallel(s2, fs, desc="Lorentz and kinematics")
    rs = " + ".join([f"({r})" for r in rs])
    # debug(f"{rs=}")

    rr = apply_den(
        rs,
        get_onshell(fds, fm) + get_mandelstamm(fds, fm),
    )
    print("len pre Optimize", len(rr))
    # TODO use #optimize from form
    rr = apply_parallel_v3(
        [rr],
        f"""
    id PREFACTOR = {"*".join([*get_color_average(fds), *get_spin_average(fds)])};
id re*i_ = 0;
id im*i_ = 1;
id re = 1;
id im = 0;
    Format {"O4" if optimize else "O0"};
    print TMP;
    .end
    """,
        desc="Optimize",
    )[0]
    # debug(f"{rr=}")
    print("len post Optimize", sum([len(relem) for relem in rr]))

    if drop_ms_prefix:
        rr = [
            (a, b.replace("ms_s", "s").replace("ms_u", "u").replace("ms_t", "t"))
            for a, b in rr
        ]

    ret = [(form.sympyfy(r[0]), form.sympyfy(r[1])) for r in rr]
    if only_result:
        if len(ret) == 1:
            return ret[0][1]
        return form.deoptimize(ret)
    return ret
