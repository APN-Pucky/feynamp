import re

import form

count = 0
dummy = 0
init = """
Symbols Pi,G,ZERO,Tr,Nc,Cf,CA,Mass,mss,mst,msu;
AutoDeclare Index Mu,Spin,Pol,Col,Glu,Propagator;
AutoDeclare Vector Mom;
Tensors f(antisymmetric),Metric(symmetric),df(symmetric),da(symmetric);
Function ProjM,ProjP,VF,xg,xgi,P,dg,dgi,xeg,xegi;
CFunctions Den,T,Denom,P,Gamma,u,v,ubar,vbar,eps,epsstar,VC,VA,GammaId;
Indices a,o,n,m,tm,tn,beta,b,m,betap,alphap,a,alpha,ind,delta,k,j,l,c,d;
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
    s = s.replace("eps_star", "epstar")
    s = s.replace("ZERO", "0")
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
