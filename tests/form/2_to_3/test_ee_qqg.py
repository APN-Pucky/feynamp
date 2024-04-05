import logging

import numpy as np
from equation_database import doi_10_1103_PhysRevD_16_3251 as ref
from feynml.interface.qgraf import style
from pyqgraf import model, qgraf

logger = logging.getLogger("feynamp")
logger.setLevel(logging.DEBUG)

import sympy
from feynmodel.interface.qgraf import feynmodel_to_qgraf
from feynmodel.interface.ufo import load_ufo_model
from pyfeyn2.feynmandiagram import FeynML
from sympy import simplify
from sympy.parsing.sympy_parser import parse_expr
from xsdata.formats.dataclass.parsers import XmlParser

from feynamp.amplitude import multiply, square_parallel
from feynamp.form import compute_squared


def test_ee_qqg():
    fm = load_ufo_model("ufo_sm")
    qfm = feynmodel_to_qgraf(fm, True, False)

    qgraf.install()
    xml_string = qgraf.run(
        "e_minus[p1], e_plus[p2]",
        "u[p3], u_bar[p4], g[p5]",
        loops=0,
        loop_momentum="l",
        model=qfm,
        style=style,
        debug=False,
    )

    parser = XmlParser()
    fml = parser.from_string(xml_string, FeynML)

    fd1, fd2 = np.array(fml.diagrams)[[f.has_pdgid(22) for f in fml.diagrams]]
    fds = [fd1, fd2]

    ret = compute_squared([fd1, fd2], fm, tag=False)

    # ret.subs("Nc","3").subs("Cf","4/3").subs("s34","-t13-t23-s35").subs("Mass_Me" , 0).simplify()

    ret = simplify(
        ret.subs("Nc", "3")
        .subs("Cf", "4/3")
        .subs("Mass_Me", 0)
        .subs("s35", "-t13-t23-s34")
        .subs("t15", "-s12-t13-t14")
        .subs("t25", "-s12-t23-t24")
        .subs("s45", "-t14-s34-t24")
    )

    # check https://www.uni-muenster.de/imperia/md/content/physik_tp/theses/klasen/neuwirth_bsc.pdf
    sret = ret.subs(
        {
            "t13": "-2*s12*x_q/4*(1-cos(theta))",
            "t23": "-2*s12*x_q/4*(1+cos(theta))",
            "t14": "-2*s12*x_qbar/4*(1-(((2-x_qbar-x_q)**2-x_qbar**2-x_q**2)/(2*x_q*x_qbar))*cos(theta)+sqrt(1-(((2-x_qbar-x_q)**2-x_qbar**2-x_q**2)/(2*x_q*x_qbar))**2)*sin(theta)*sin(eta))",
            "t24": "-2*s12*x_qbar/4*(1+(((2-x_qbar-x_q)**2-x_qbar**2-x_q**2)/(2*x_q*x_qbar))*cos(theta)-sqrt(1-(((2-x_qbar-x_q)**2-x_qbar**2-x_q**2)/(2*x_q*x_qbar))**2)*sin(theta)*sin(eta))",
            "s34": "2*s12/4*x_q*x_qbar * (1-(((2-x_qbar-x_q)**2-x_qbar**2-x_q**2)/(2*x_q*x_qbar)))",
        }
    ).simplify()

    int2 = sympy.integrate(
        sret * sympy.parse_expr("sin(theta)"), ("theta", 0, sympy.pi)
    )
    int1 = sympy.integrate(int2, ("eta", 0, 2 * sympy.pi))
    int3 = sympy.integrate(int1, ("phi", 0, 2 * sympy.pi))
    ssret = int3.simplify()

    # integration prefactors
    fret = ssret / ((2 * sympy.pi) ** 5 * 64)

    assert fret.equals(
        (ref.equation_2_9.rhs * ref.equation_2_8.rhs).subs(
            {
                "x_1": "x_q",
                "x_2": "x_qbar",
                "alpha_C": "G**2/4/pi",
                "alpha": "ee**2/4/pi",
                "Q": "s12**(1/2)",
                "sum_e_q_squared": "(2/3)**2",
            }
        )
    )
