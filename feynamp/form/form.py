import os
import re

import form

count = 0
dummy = 0
init = """
Symbols Pi,G,ZERO,Tr,Nc,Cf,CA,MC,ee;
AutoDeclare Index Mu,Spin,Pol,Col,Glu,Propagator;
AutoDeclare Symbol Mass,fd,mss,mst,msu;
AutoDeclare Vector Mom;
Tensors f(antisymmetric),Metric(symmetric),df(symmetric),da(symmetric),Identity(symmetric);
Function ProjM,ProjP,VF,xg,xgi,P,dg,dgi,xeg,xegi;
CFunctions Den,T,Denom,P,Gamma,u,v,ubar,vbar,eps,epsstar,VC,VA,GammaId, GammaCollect, GammaIdCollect;
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


def run_parallel(init, cmds, variables, show=False, keep_form_file=False, threads=None):
    global count
    count = count + 1
    rets = []
    if threads is None:
        threads = os.cpu_count()
        with form.open(keep_log=1000, args=["tform", f"-w{threads}"]) as f:
            txt = "" + init
            for i, s in enumerate(variables):
                txt += f"Local TMP{i} = {s};\n"
            txt += cmds
            for i, s in enumerate(variables):
                txt += f"print TMP{i};.sort;"
            f.write(txt)
            if keep_form_file:
                with open("form" + str(count) + ".frm", "w") as frm:
                    frm.write(txt)
            for i, s in enumerate(variables):
                rets.append(f.read(f"TMP{i}"))
            # What is this ?
            # r = re.sub(r"\+factor_\^?[0-9]*", r"", r).strip("*")
            if show:
                for r in rets:
                    print(r + "\n")
            assert len(rets) == len(variables)
            return rets


def run(s, show=False, keep_form_file=False):
    global count
    count = count + 1
    with form.open(keep_log=1000, args=["tform", "-w16"]) as f:
        local = s.split("Local")[1].split("=")[0].strip()
        txt = s + "print " + local + ";.sort;"
        f.write(txt)
        # frm.write(txt)
        r = f.read("" + local)
        r = re.sub(r"\+factor_\^?[0-9]*", r"", r).strip("*")
        if show:
            print(r + "\n")
        # delete frm file
        if keep_form_file:
            with open("form" + str(count) + ".frm", "w") as frm:
                frm.write(txt)
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
