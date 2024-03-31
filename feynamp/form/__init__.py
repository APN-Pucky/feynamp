import os
import re
from typing import List

import form
import sympy
from feynml.feynmandiagram import FeynmanDiagram
from feynmodel.feyn_model import FeynModel

import feynamp
from feynamp.log import debug

count = 0
dummy = 0
init = """
Symbols Pi,G,ZERO,Tr,Nc,Cf,CA,MC,ee;
AutoDeclare Index Mu,Spin,Pol,Col,Glu,Propagator;
AutoDeclare Symbol Mass,fd,mss,mst,msu;
AutoDeclare Vector Mom;
Tensors f(antisymmetric),Metric(symmetric),df(symmetric),da(symmetric),Identity(symmetric);
Function ProjM,ProjP,VF,xg,xgi,P,dg,dgi,xeg,xegi;
CFunctions Den,T,Denom,P,Gamma,u,v,ubar,vbar,eps,epsstar,VC,VA,GammaId;
Indices a,o,n,m,tm,tn,beta,b,m,betap,alphap,a,alpha,ind,delta,k,j,l,c,d,e;
"""


def get_dummy_index():
    global dummy
    dummy = dummy + 1
    return f"N{dummy}_?"


def string_to_form(s):
    s = s.replace("complex(0,1)", "i_")  # form uses i_ for imaginary unit
    s = s.replace("Gamma_Id", "GammaId")
    s = s.replace("u_bar", "ubar")
    s = s.replace("v_bar", "vbar")
    s = s.replace("eps_star", "epsstar")
    s = s.replace("Identity", "df")  # TODO check if this holds or also happens for anti
    s = s.replace("ZERO", "0")
    s = s.replace(".*", "*")  # handle decimals
    s = s.replace(".)", ")")  # handle decimals
    return s


def run(s, show=False, keep_form_file=True):
    global count
    count = count + 1
    with open("form" + str(count) + ".frm", "w") as frm:
        with form.open(keep_log=1000) as f:
            l = s.split("Local")[1].split("=")[0].strip()
            txt = s + "print " + l + ";.sort;"
            f.write(txt)
            frm.write(txt)
            r = f.read("" + l)
            r = re.sub(r"\+factor_\^?[0-9]*", r"", r).strip("*")
            if show:
                print(r + "\n")
            # delete frm file
            if not keep_form_file:
                os.remove("form" + str(count) + ".frm")
            return r


def sympyfy(string_expr):
    from sympy import simplify
    from sympy.parsing.sympy_parser import parse_expr

    ret = simplify(
        parse_expr(
            string_expr.replace("Mom_", "")
            .replace(".", "_")
            .replace("^", "**")
            .replace("mss", "s")
            .replace("msu", "u")
            .replace("mst", "t")
        )
    )
    return simplify(ret.subs("Nc", "3").subs("Cf", "4/3"))


# TODO compute squared  functino which coutns legs!!"!!!" and picks right mandelstamm,s


def compute_squared(fds: List[FeynmanDiagram], fm: FeynModel):
    dims = fds[0].get_externals_size()
    for fd in fds:
        assert (
            dims == fd.get_externals_size()
        ), "All FeynmanDiagrams must have the same external legs"
    s2 = feynamp.amplitude.square(fds, fm, tag=False)
    debug(f"{s2=}")
    fs = ""
    fs += feynamp.form.lorentz.get_gammas()
    fs += feynamp.form.color.get_color()
    fs += feynamp.form.momentum.get_kinematics()
    fs += feynamp.form.momentum.get_onshell(fds, fm)
    fs += feynamp.form.momentum.get_mandelstamm(fds, fm)

    rs = feynamp.form.momentum.apply(s2, fs)
    debug(f"{rs=}")

    rr = feynamp.form.momentum.apply_den(
        rs,
        feynamp.form.momentum.get_onshell(fds, fm)
        + feynamp.form.momentum.get_mandelstamm(fds, fm),
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
        ).subs({"Nc": "3", "Cf": "4/3"})
        * sympy.parse_expr(
            "*".join([*feynamp.get_color_average(fds), *feynamp.get_spin_average(fds)])
        )
    )
    return ret
