import os
import re
import subprocess
import tempfile
from typing import List

import form
import tqdm
from pqdm.threads import pqdm

from feynamp.leg import get_leg_momentum

count = 0
dummy = 0
# TODO auto generate symbols
init = """
Symbols Pi,G,ZERO,Tr,Nc,Cf,CA,MC,ee,realpart,PREFACTOR;
AutoDeclare Index Mu,Spin,Pol,Propagator;
AutoDeclare Symbol Mass,fd;
* Mandelstamm
AutoDeclare Symbol ms;
* Momentum
AutoDeclare Vector Mom;
Tensors colorcorrelation,spincorrelation;
Index scMuMu,scMuNu;
Tensors Metric(symmetric),df(symmetric),da(symmetric),Identity(symmetric);
Function ProjM,ProjP,VF,xg,xgi,P,dg,dgi,xeg,xegi;
CFunctions Den,Denom,P,Gamma,u,v,ubar,vbar,eps,epsstar,VC,VA,VPol,GammaId, GammaCollect, GammaIdCollect;
Indices a,o,n,m,tm,tn,beta,b,m,betap,alphap,a,alpha,ind,delta,k,j,l,c,d,e;
"""


def get_dummy_index(underscore=True, questionmark=True):
    global dummy
    dummy = dummy + 1
    return f"N{dummy}" + ("_" if underscore else "") + ("?" if questionmark else "")


def string_to_form(s):
    try:
        s = re.sub(r"complex\((.*?),(.*?)\)", r"(\1+i_*(\2))", s)
        # s = s.replace("complex(0,1)", "i_")  # form uses i_ for imaginary unit
        s = s.replace("Gamma_Id", "GammaId")
        s = s.replace("u_bar", "ubar")
        s = s.replace("v_bar", "vbar")
        s = s.replace("eps_star", "epsstar")
        s = s.replace(
            "Identity", "df"
        )  # TODO check if this holds or also happens for anti
        s = s.replace("ZERO", "0")
        s = s.replace(".*", "*")  # handle decimals
        s = s.replace(".)", ")")  # handle decimals
    except Exception as e:
        print("Error in string_to_form", e)
        print(s)
        raise e
    return s


def apply_parallel_v3(amps: List[str], operations: str, desc=None):
    return run_parallel_v3(
        init, operations, [string_to_form(a) for a in amps], desc=desc
    )


def run_parallel(*args, **kwargs):
    # return run_parallel_v1(*args, **kwargs)
    return run_parallel_v2(*args, **kwargs)


def run_parallel_v3(
    init, cmds, variables, show=False, keep_form_file=True, threads=None, desc=None
):
    global count
    count = count + 1
    rets = []
    if threads is None:
        threads = os.cpu_count()
    ret = pqdm(
        [{"tuples": [("TMP", var)], "init": init, "code": cmds} for var in variables],
        run_bare_multi,
        n_jobs=threads,
        desc=desc,
        argument_type="kwargs",
    )
    # check if ret is error
    if isinstance(ret, Exception):
        raise ret
    if isinstance(ret, List):
        for r in ret:
            if isinstance(r, Exception):
                raise r
    return ret


def run_parallel_v2(
    init, cmds, variables, show=False, keep_form_file=True, threads=None, desc=None
):
    global count
    count = count + 1
    rets = []
    if threads is None:
        threads = os.cpu_count()
    return pqdm(
        ["" + init + f"Local TMP = {var};" + cmds for var in variables],
        run_bare,
        n_jobs=threads,
        desc=desc,
    )


def run_parallel_v1(
    init, cmds, variables, show=False, keep_form_file=True, threads=None
):
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
            # Not sure why sort is needed, but it is needed
            txt += f"print TMP{i};.sort;"
        if keep_form_file:
            with open("form" + str(count) + ".frm", "w") as frm:
                frm.write(txt)
        f.write(txt)
        for i, s in enumerate(variables):
            rets.append(f.read(f"TMP{i}"))
        # What is this ?
        # r = re.sub(r"\+factor_\^?[0-9]*", r"", r).strip("*")
        if show:
            for r in rets:
                print(r + "\n")
        assert len(rets) == len(variables)
        return rets


def run_bare(s, show=False, keep_form_file=True):
    init = s.split("Local")[0]
    local = s.split("Local")[1].split("=")[0].strip()
    eq = s.split("Local")[1].split("=")[1].strip().split(";")[0]
    code = ";".join("=".join(s.split("Local")[1].split("=")[1:]).split(";")[1:])
    ret = run_bare_multi(
        [(local, eq)],
        init=init,
        code=code,
        show=show,
        keep_form_file=keep_form_file,
        end="print " + local + ";\n.sort\n.end",
    )
    if len(ret) != 1:
        raise ValueError(f"Error1 in form output. found #{len(ret)} ")
    if len(ret[0]) != 2:
        raise ValueError(f"Error2 in form output. found #{len(ret[0])} ")
    return ret[0][1]


def run_bare_multi(tuples, init, code, show=False, keep_form_file=True, end=".end"):
    """Run it just as a subprocess"""
    # make temporary file
    with tempfile.NamedTemporaryFile(
        "w", suffix=".frm", delete=not keep_form_file
    ) as f:
        txt = "Off Statistics;\n"
        txt += init
        for s in tuples:
            txt += f"Local {s[0]} = {s[1]};\n"
        txt += code
        txt += end
        f.write(txt)
        # flush it
        f.flush()
        # run form on file and capture output
        out = subprocess.check_output(["form", "-q", f.name]).decode()
        sout = out.replace("\n", "").replace(" ", "")
        res = re.findall(r"(.*?)=(.*?);", sout, re.DOTALL)
        # print(res)
        return res


def run(s, show=False, keep_form_file=True, threads=1):
    global count
    count = count + 1
    if threads is None:
        threads = os.cpu_count()
    with form.open(keep_log=1000, args=["tform", f"-w{threads}"]) as f:
        local = s.split("Local")[1].split("=")[0].strip()
        # no clue why sort
        txt = s + "print " + local + ";.sort;"
        if keep_form_file:
            with open("form" + str(count) + ".frm", "w") as frm:
                frm.write(txt)
        f.write(txt)
        # frm.write(txt)
        r = f.read("" + local)
        r = re.sub(r"\+factor_\^?[0-9]*", r"", r).strip("*")
        if show:
            print(r + "\n")
        return r


def deoptimize(ret):
    import sympy

    symbols = []
    for r in ret:
        symbols += list(r[1].free_symbols)
    sss = set(symbols)
    ss = ",".join([str(s) for s in sss])
    exec(f'{ss} = sympy.symbols("{ss}")')
    for r in tqdm.tqdm(ret, desc="Deoptimizing"):
        exec(f"{r[0]}={r[1]}")
    return eval(str(ret[-1][0]))


def sympyfy(string_expr):
    from sympy import simplify
    from sympy.parsing.sympy_parser import parse_expr

    # ret = simplify(
    ret = parse_expr(
        string_expr
        # .replace("Mom_", "")
        .replace(".", "_").replace("^", "**")
        # .replace("ms_s", "s")
        # .replace("ms_u", "u")
        # .replace("ms_t", "t")
        ,
        evaluate=False,
    )
    # )
    return ret
    # return simplify(ret.subs("Nc", "3").subs("Cf", "4/3"))


def sympy_to_form_string(sympy_expr):
    return str(sympy_expr).replace("**", "^")
