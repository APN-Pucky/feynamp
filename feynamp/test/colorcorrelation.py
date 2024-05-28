from feynamp import sympy
from feynamp.leg import color_vector_to_casimir, get_color_vector, get_leg_momentum


def assert_colorcorrelation(sympy_expr, fds, model):
    """
    sympy_expr: sympy expression of cc/born
    """
    fd = fds[0]
    legs = fd.legs
    for j in range(len(legs)):
        vec = get_color_vector(fd, legs[j], model)
        if vec is not None:
            casimir = color_vector_to_casimir(vec)
            sum = sympy_expr
            for i in range(len(legs)):
                momi = get_leg_momentum(legs[i])
                momj = get_leg_momentum(legs[j])
                # colorcorrelation is symmetric, and we only have sorted vertices
                sum = sum.replace(
                    sympy.parse_expr(f"colorcorrelation({momi},{momj})"),
                    1,
                )
                sum = sum.replace(
                    sympy.parse_expr(f"colorcorrelation({momj},{momi})"),
                    1,
                )
            # TODO can be optimized
            for k in range(len(legs)):
                for l in range(len(legs)):
                    momi = get_leg_momentum(legs[k])
                    momj = get_leg_momentum(legs[l])
                    sum = sum.replace(
                        sympy.parse_expr(
                            f"colorcorrelation({momi},{momj})".replace("Mom_", "")
                        ),
                        0,
                    )
                    sum = sum.replace(
                        sympy.parse_expr(
                            f"colorcorrelation({momj},{momi})".replace("Mom_", "")
                        ),
                        0,
                    )
            # replace all remaining colorcorrelation(mom,mom) with 0
            assert sum.equals(sympy.parse_expr(casimir))
