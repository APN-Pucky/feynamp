from feynamp.form import *

color = """
**********************************************************
*                  COLOUR STRUCTURE SIMPLIFY             *
**********************************************************
    
repeat;
* remove df(k,j)
   id df(k?,l?)*df(l?,j?)=df(k,j);
   id T(a?,k?,l?)*df(k?,j?)=T(a,j,l);id T(a?,k?,l?)*df(l?,j?)=T(a,k,j);
* remove da(a,b)
   id da(a?,b?)*da(b?,c?)=da(a,c);
   id T(a?,k?,l?)*da(a?,b?)=T(b,k,l);
   id f(a?,b?,c?)*da(a?,d?)=f(d,b,c);
* length-three objects simplify:
   id T(b?,k?,j?)*T(a?,j?,c?)*T(b?,c?,l?)=(Cf-Nc*Tr)*T(a,k,l);
   id T(b?,j?,l?)*T(c?,l?,k?)*f(a?,b?,c?)=i_*Nc*Tr*T(a,j,k);
* length-two objects that give out df(k,j)
   id T(a?,c?,j?)*T(a?,k?,l?)=-1/Nc*df(c, k)*df(j, l)/2 + df(c, l)*df(j, k)/2;
   id T(a?,k?,l?)*T(a?,l?,j?)=Cf*df(k,j);
* length-two objects that give out da(a,b)
   id T(a?,k?,l?)*T(b?,l?,k?)=Tr*da(a,b);
   id f(a?,b?,c?)*f(d?,b?,c?)=Nc*da(a,d); 
* simplify traces
   id T(b?,k?,k?)=0;
   id da(a?,a?)=Nc*Cf/Tr;
   id df(a?,a?)=Nc;
* simplify combination of factors
   id Nc^-2=2-Nc^2+Cf^2*Tr^-2;
   id Nc^2=1+Nc*Cf/Tr;
   id Tr=1/2;
   id Tr^-1=2;
* double f(a,b,c) simplify
*  id f(a?,b?,e?)*f(c?,d?,e?)=-2 * { [T(a), T(b)] [T(c), T(d)]};
   id f(a?,b?,e?)*f(c?,d?,e?)=-2 * (T(a,e,N1_?)*T(b,N1_?,N2_?) - T(b,e,N1_?)*T(a,N1_?,N2_?))*(T(c,N2_?,N3_?)*T(d,N3_?,e)-T(d,N2_?,N3_?)*T(c,N3_?,e));
endrepeat;
"""

color_sum = """
**********************************************************
*                  COLOUR SUM SIMPLIFY                   *
**********************************************************
repeat;
  id VA(Glua?,Momb?)*VA(Gluc?,Momb?) = da(Glua,Gluc);
  id VC(Cola?,Momb?)*VC(Colc?,Momb?) = df(Cola,Colc);
endrepeat;
"""


def get_color():
    return color_sum + color


def apply_color(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + get_color())


def apply_color_sum(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + color_sum)


def apply_color_simplify(string_expr):
    s = string_to_form(string_expr)
    return run(init + f"Local TMP = {s};" + color)
