import re

import form

count = 0
init = """
Symbols Pi,G,ZERO,Tr,Nc,Cf,CA;
AutoDeclare Index Mu,Mom,Spin,Pol,Col,Glu,Propagator;
Tensors f(antisymmetric),Metric(symmetric),df(symmetric),da(symmetric);
Function ProjM,ProjP,VF,xg,xgi,P,dg,dgi,xeg,xegi;
CFunctions T,Denom,P,Gamma,u,v,ubar,vbar,eps,epsstar,VC,VA;
Indices a,o,n,m,tm,tn,beta,b,betap,alphap,a,alpha,ind,delta,k,j,l,c,d;
"""


def string_to_form(s):
    s = s.replace("complex(0,1)", "i_")  # form uses i_ for imaginary unit
    s = s.replace("u_bar", "ubar")
    s = s.replace("v_bar", "vbar")
    s = s.replace("eps_star", "epstar")
    return s


def run(s, show=False):
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
            return r
