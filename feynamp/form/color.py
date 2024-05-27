from typing import List

from feynamp.form.form import get_dummy_index, init, run, run_parallel, string_to_form
from feynamp.leg import (
    color_vector_to_id,
    color_vector_to_index,
    color_vector_to_operator,
    get_color_vector,
    get_leg_momentum,
    is_swapped_color_vector,
)
from feynamp.log import debug

from .colorh import colorh

colorh_init = f"""
{colorh}
AutoDeclare Index Color=NR;
AutoDeclare Index Glu=NA;
"""

color_init = """
AutoDeclare Index Color;
AutoDeclare Index Glu;
Tensors f(antisymmetric), colorcorrelation;
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
* https://www.nikhef.nl/~form/maindir/packages/color/color.pdf
* Appendix D
id d33(cOlpR1,cOlpR2) = I2R^3/(2*Nc)*(Nc^2-1)*(Nc^2-4);
id d44(cOlpR1,cOlpR2) = I2R^4/(6*Nc^2)*(Nc^2-1)*(Nc^4-6*Nc^2+18);
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


def get_color(fds=None, legs=None, model=None, colorcorrelation=False):
    assert not colorcorrelation or (
        fds is not None and legs is not None and model is not None
    )
    # return get_color_v1()
    return get_color_v3(fds, legs, model, colorcorrelation=colorcorrelation)


def get_full_color_correlation_matrix(fds, legs, model):

    left = ""
    right = ""
    vec = []
    ops = []
    ind = []
    ind1 = []
    ind2 = []
    mom = []
    swap = []
    ids = []
    for i in range(len(legs)):
        swap += [is_swapped_color_vector(fds[0], legs[i], model)]
        vec += [get_color_vector(fds[0], legs[i], model)]
        ids += [color_vector_to_id(vec[i])]
        mom += [get_leg_momentum(legs[i])]
        ops += [color_vector_to_operator(vec[i])]
        ind += [color_vector_to_index(vec[i])]
        ind1 += [str(ind[i]) + legs[i].id]
        ind2 += [str(ind[i]) + get_dummy_index(underscore=False, questionmark=False)]
    dummy = "Glu" + get_dummy_index(underscore=False, questionmark=False)
    for i in range(len(legs)):
        if vec[i] is not None:
            i1 = ind1[i] + " "
            i2 = ind2[i] + "?"
            if swap[i]:
                i1, i2 = i2, i1
            left += f"{vec[i]}({i1},{mom[i]})*{vec[i]}({i2},{mom[i]})*"
            i1, i2 = i1[:-1], i2[:-1]
            for j in range(i + 1, len(legs)):
                if vec[j] is not None:

                    j1 = ind1[j]
                    j2 = ind2[j]
                    if swap[j]:
                        j1, j2 = j2, j1
                    deltas = "*"
                    for k in range(len(legs)):
                        if vec[k] is not None and i != k and j != k:
                            deltas += f"{ids[k]}({ind1[k]},{ind2[k]})*"
                    right += f"\ncolorcorrelation({mom[i]},{mom[j]})*{ops[i]}({i1},{i2},{dummy})*{ops[j]}({j1},{j2},{dummy}){deltas[:-1]}+"
                    # right += f"\n{ops[i]}({i1},{i2},{dummy})*{ops[j]}({j1},{j2},{dummy}){deltas[:-1]}+"
                    # right += f"\ncolorcorrelation({mom[i]},{mom[j]})*d_({i1},{i2})*d_({j1},{j2}){deltas[:-1]}+"
    ret = f"""
    id {left[:-1]} = {right[:-1]};
    """
    print(ret)
    return ret


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


def get_color_v3(fds, legs, model, colorcorrelation=False):
    ret = ""
    if colorcorrelation:
        ret += get_full_color_correlation_matrix(fds, legs, model)
    else:
        ret += color_sum
    return ret + colorh_ids + color + colorh_ids + rep(color_simplify)


def get_color_v2():
    return color_sum + colorh_ids + color + colorh_ids + rep(color_simplify)


def get_color_v1():
    return color_sum + color_ids + old_color


def get_color_ids():
    return color_sum + colorh_ids


def apply_color_ids(string_expr):
    s = string_to_form(string_expr)
    return run(init + colorh_init + f"Local TMP = {s};" + get_color_ids())


def apply_color_parallel(string_exprs: List[str], **kwargs):
    return run_parallel(
        init + colorh_init,
        get_color(**kwargs),
        [string_to_form(a) for a in string_exprs],
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
