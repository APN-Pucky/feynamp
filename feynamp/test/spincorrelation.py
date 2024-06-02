import sympy

from feynamp.form.form import sympy_to_form_string, sympyfy
from feynamp.form.momentum import apply, get_kinematics, get_mandelstamm, get_onshell
from feynamp.leg import (
    color_vector_to_casimir,
    get_color_vector,
    get_leg_momentum,
    is_vector,
)


def assert_spincorrelation(sympy_expr, fds, model):
    """
    sympy_expr: sympy expression of sc/born

    Checks https://arxiv.org/pdf/1002.2581 Eq. 2.10
    """
    fd = fds[0]
    legs = fd.legs
    for j, leg in enumerate(legs):
        if is_vector(fd, leg, model):
            mom = get_leg_momentum(leg)
            sum = sympy_expr
            sum = sum.replace(
                sympy.parse_expr(f"spincorrelation({mom},scMuMu,scMuNu)"), 1
            )
            for i, leg in enumerate(legs):
                mom = get_leg_momentum(leg)
                sum = sum.replace(
                    sympy.parse_expr(f"spincorrelation({mom},scMuMu,scMuNu)"), 0
                )
            sum = sum * sympy.parse_expr("d_(MuMu,MuNu)")
            st = sympy_to_form_string(sum)
            fs = ""
            fs += get_kinematics()
            fs += get_onshell(fds, model)
            fs += get_mandelstamm(fds, model)

            rr = apply(st, fs)
            sr = sympyfy(rr).simplify()
            print(sr)
            assert sr == sympy.parse_expr("-1")
