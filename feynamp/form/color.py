from typing import List

from feynamp.form.form import get_dummy_index, init, run, run_parallel, string_to_form

from .colorh import colorh

colorh_init = f"""
{colorh}
AutoDeclare Index Color=NR;
AutoDeclare Index Glu=NA;
"""

color_init = """
AutoDeclare Index Color;
AutoDeclare Index Glu;
Tensors f(antisymmetric);
CFunctions T;
Symbols NA,I2R;
"""

color_ids = """
repeat;
* remove df(k,j)
   id df(k?,l?)*df(l?,j?)=df(k,j);
   id T(a?,k?,l?)*df(k?,j?)=T(a,j,l);
   id T(a?,k?,l?)*df(l?,j?)=T(a,k,j);
* remove da(a,b)
   id da(a?,b?)*da(b?,c?)=da(a,c);
   id T(a?,k?,l?)*da(a?,b?)=T(b,k,l);
   id f(a?,b?,c?)*da(a?,d?)=f(d,b,c);
* simplify traces
   id T(b?,k?,k?)=0;
   id da(a?,a?)=Nc*Cf/Tr;
   id df(a?,a?)=Nc;
endrepeat;
"""

colorh_ids = """
repeat;
* remove df(k,j)
   id df(k?,l?)*df(l?,j?)=df(k,j);
   id T(k?,l?,a?)*df(k?,j?)=T(j,l,a);
   id T(k?,l?,a?)*df(l?,j?)=T(k,j,a);
* remove da(a,b)
   id da(a?,b?)*da(b?,c?)=da(a,c);
   id T(k?,l?,a?)*da(a?,b?)=T(k,l,b);
   id f(b?,c?,a?)*da(a?,d?)=f(b,c,d);
* simplify traces
   id T(k?,k?,b?)=0;
   id da(a?,a?)=Nc*Cf/Tr;
   id df(a?,a?)=Nc;
endrepeat;
"""

color_simplify = """
* simplify combination of factors
*   id Nc^-2=2-Nc^2+Cf^2*Tr^-2;
*   id Nc^2=1+Nc*Cf/Tr;
   id NA=Nc*Cf/Tr;
   id cA=Nc;
   id cR=Cf;
   id I2R=1/2;
   id Tr=1/2;
   id Tr^-1=2;
"""

old_color = f"""
**********************************************************
*                  COLOUR STRUCTURE SIMPLIFY             *
**********************************************************
    
repeat;
{color_ids}
* length-three objects simplify:
   id T(b?,k?,j?)*T(a?,j?,c?)*T(b?,c?,l?)=(-Tr/Nc*T(a,k,l));
   id T(b?,j?,l?)*T(c?,l?,k?)*f(a?,b?,c?)=(i_*Nc*Tr*T(a,j,k));
* length-two objects that give out df(k,j)
   id T(a?,c?,j?)*T(a?,k?,l?)=(-1/Nc*df(c, j)*df(k, l)*Tr + df(c, l)*df(j, k)*Tr);
   id T(a?,k?,l?)*T(a?,l?,j?)=(Cf*df(k,j));
* length-two objects that give out da(a,b)
   id T(a?,k?,l?)*T(b?,l?,k?)=(Tr*da(a,b));
   id f(a?,b?,c?)*f(d?,b?,c?)=Nc*da(a,d); 
{color_simplify}
* double f(a,b,c) simplify
*  id f(a?,b?,e?)*f(c?,d?,e?)=-2 * {{ [T(a), T(b)] [T(c), T(d)]}};
   id f(a?,b?,e?)*f(c?,d?,e?)=-2 * (T(a,e,N1_?)*T(b,N1_?,N2_?) - T(b,e,N1_?)*T(a,N1_?,N2_?))*(T(c,N2_?,N3_?)*T(d,N3_?,e)-T(d,N2_?,N3_?)*T(c,N3_?,e));
endrepeat;
"""


color = """
#call docolor
"""
# TODO do the color stuff manually (cf. MG) since the simplifications here are very expensive

color_sum = """
**********************************************************
*                  COLOUR SUM SIMPLIFY                   *
**********************************************************
repeat;
  id VA(Glua?,Momb?)*VA(Gluc?,Momb?) = da(Glua,Gluc);
  id VC(Colora?,Momb?)*VC(Colorc?,Momb?) = df(Colora,Colorc);
endrepeat;
"""


def rep(s: str):
    return f"""
repeat;
   {s}
endrepeat;
"""


repeat = rep


def get_color(mom1=None, mom2=None):
    # return get_color_v1()
    return get_color_v3(mom1=mom1, mom2=mom2)


def get_color_sum_v1(mom1=None, mom2=None):
    ret = ""
    if mom1 is not None and mom2 is not None:
        # new index that is summed over
        dummy = "Glu" + get_dummy_index()
        ret += repeat(
            f"""
         id VA(Glua?,{mom1})*VA(Gluc?,{mom1}) = f(Glua,Gluc,{dummy});
         id VA(Glua?,{mom2})*VA(Gluc?,{mom2}) = f(Glua,Gluc,{dummy});
         id VC(Colora?,{mom1})*VC(Colorc?,{mom1}) = T(Colora,Colorc,{dummy});
         id VC(Colora?,{mom2})*VC(Colorc?,{mom2}) = T(Colora,Colorc,{dummy});
"""
        )
    return ret + color_sum


def get_color_v3(mom1=None, mom2=None):
    return (
        get_color_sum_v1(mom1=mom1, mom2=mom2)
        + colorh_ids
        + color
        + colorh_ids
        + rep(color_simplify)
    )


def get_color_v2():
    return color_sum + colorh_ids + color + colorh_ids + rep(color_simplify)


def get_color_v1():
    return color_sum + color_ids + old_color


def get_color_ids():
    return color_sum + colorh_ids


def apply_color_ids(string_expr):
    s = string_to_form(string_expr)
    return run(init + colorh_init + f"Local TMP = {s};" + get_color_ids())


def apply_color_parallel(string_exprs: List[str]):
    return run_parallel(
        init + colorh_init, get_color(), [string_to_form(a) for a in string_exprs]
    )


def apply_color(string_expr):
    s = string_to_form(string_expr)
    return run(init + colorh_init + f"Local TMP = {s};" + get_color())


def apply_color_sum(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + color_sum)


def apply_color_simplify(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + color)
