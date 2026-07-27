#!/usr/bin/env python3
"""Recover and certify the universal conjugate-cubic norm identity.

The computation is over K = Q(lam,e2,e3,u) on the e1=1, v=1 chart.  The
sparse b2 and b1 candidates were discovered by exact interpolation, but no
specialization is used to prove the norm identity: all seven coefficients are
checked in K.

The monic gcd assertion for the two residual b1 equations is certified as
follows.  The universal identity proves that x-b1 divides both equations.  At
the exact point (lam,e2,e3,u)=(3,2,5,1), their degrees remain 3 and 4 and their
gcd has degree one.  If the generic gcd had degree at least two, the relevant
subresultants would vanish identically and hence at this good specialization.
Thus the generic monic gcd is exactly x-b1.

The saved formulas are compared read-only by default.  The explicit
--write-formulas option replaces them after the identity has been verified.
"""

from __future__ import annotations

import argparse
import pathlib
import re

import sympy as sp


HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "generated" / "g2_mu6_universal_cover_e1one_vone.sing"
FORMULAS = HERE / "generated" / "g2_mu6_e1one_vone_norm_formulas.txt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-formulas",
        action="store_true",
        help="replace the saved formula file instead of verifying it read-only",
    )
    return parser.parse_args()


def read_assignment(name: str) -> str:
    text = SOURCE.read_text(encoding="ascii")
    match = re.search(rf"^poly {re.escape(name)}=(.*);$", text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"missing {name} assignment in {SOURCE}")
    return match.group(1).replace("^", "**")


def main(*, write_formulas: bool = False) -> None:
    lam, e2, e3, u, z, x = sp.symbols("lam e2 e3 u z x")
    parameters = (lam, e2, e3, u)
    local_dict = {str(symbol): symbol for symbol in (*parameters, z)}
    coefficient_field = sp.QQ.frac_field(*parameters)

    curve = sp.sympify(read_assignment("curve"), locals=local_dict)
    four_s = sp.sympify(read_assignment("fourS"), locals=local_dict)
    curve_poly = sp.Poly(curve, z, domain=coefficient_field)
    if curve_poly.degree() != 6:
        raise RuntimeError(f"expected a sextic residual equation, got degree {curve_poly.degree()}")
    leading = coefficient_field.from_sympy(curve_poly.LC())
    coefficients = [
        coefficient_field.from_sympy(curve_poly.nth(degree)) / leading
        for degree in range(7)
    ]
    if coefficients[6] != coefficient_field.one:
        raise RuntimeError("the normalized residual equation is not monic")
    print("MONIC_CURVE_DEGREE", curve_poly.degree())
    k = coefficient_field.from_sympy(3 * four_s)
    leading_numerator = sp.fraction(
        sp.cancel(coefficient_field.to_sympy(leading))
    )[0]

    b2_numerator = e2 * lam * (
        2 * e2**3 * lam * u
        + e2**3 * lam
        - 6 * e2 * e3 * lam * u
        + 6 * e2 * e3 * lam
        - 3 * e3 * lam
        + 3 * e3
    )
    b1_numerator = -lam * (
        e2**3 * lam * u
        + 2 * e2**3 * lam
        + 3 * e2**2 * e3 * lam * u**2
        - 6 * e2**2 * e3 * lam * u
        + 3 * e2**2 * e3 * lam
        - 3 * e2 * e3 * lam * u
        + 3 * e2 * e3 * lam
        - 3 * e3 * lam
        + 3 * e3
    )
    b2 = coefficient_field.from_sympy(b2_numerator / (2 * leading_numerator))
    b1 = coefficient_field.from_sympy(b1_numerator / (2 * leading_numerator))
    if b2 == coefficient_field.zero:
        raise RuntimeError("the proposed quadratic B has zero leading coefficient")

    a2 = coefficients[5] / 2
    a1 = (coefficients[4] - a2**2 - k * b2**2) / 2
    a0 = coefficients[3] / 2 - a2 * a1 - k * b2 * b1
    b0 = (
        coefficients[2] - a1**2 - 2 * a2 * a0 - k * b1**2
    ) / (2 * k * b2)

    a_coefficients = (a0, a1, a2, coefficient_field.one)
    b_coefficients = (b0, b1, b2)
    print("A_DEGREE", len(a_coefficients) - 1)
    print("B_DEGREE", len(b_coefficients) - 1)
    print("B_LEADING_COEFFICIENT_NONZERO", b2 != coefficient_field.zero)
    product_coefficients = [coefficient_field.zero for _ in range(7)]
    for left_degree, left in enumerate(a_coefficients):
        for right_degree, right in enumerate(a_coefficients):
            product_coefficients[left_degree + right_degree] += left * right
    for left_degree, left in enumerate(b_coefficients):
        for right_degree, right in enumerate(b_coefficients):
            product_coefficients[left_degree + right_degree] += k * left * right

    remainders = [
        product_coefficients[degree] - coefficients[degree]
        for degree in range(7)
    ]
    print("FIELD_COEFFICIENT_REMAINDERS", [remainder == 0 for remainder in remainders])
    if any(remainders):
        raise RuntimeError("the universal norm identity failed")

    # A single exact good specialization bounds the generic gcd degree above.
    point = {lam: sp.Integer(3), e2: sp.Integer(2), e3: sp.Integer(5), u: sp.Integer(1)}
    specialized_curve = sp.Poly(curve.subs(point), z, domain=sp.QQ).monic()
    specialized_coefficients = [specialized_curve.nth(degree) for degree in range(7)]
    specialized_k = sp.Rational((3 * four_s).subs(point))
    specialized_b2 = sp.Rational((b2_numerator / (2 * leading_numerator)).subs(point))
    specialized_b1 = sp.Rational((b1_numerator / (2 * leading_numerator)).subs(point))
    specialized_a2 = specialized_coefficients[5] / 2
    specialized_a1 = (
        specialized_coefficients[4]
        - specialized_a2**2
        - specialized_k * specialized_b2**2
    ) / 2
    p = specialized_coefficients[3] / 2 - specialized_a2 * specialized_a1
    q = specialized_k * specialized_b2
    r = specialized_coefficients[2] - specialized_a1**2 - 2 * specialized_a2 * p
    d0 = r / (2 * q)
    d2 = -1 / (2 * specialized_b2)
    equation_one = sp.Poly(
        2 * specialized_a1 * p
        - specialized_coefficients[1]
        + (-2 * specialized_a1 * q + r / specialized_b2) * x
        + 2 * specialized_k * specialized_a2 * x**2
        - specialized_k / specialized_b2 * x**3,
        x,
        domain=sp.QQ,
    )
    equation_zero = sp.Poly(
        p**2
        + specialized_k * d0**2
        - specialized_coefficients[0]
        + (-2 * p * q + 2 * specialized_k * d0 * specialized_a2) * x
        + (
            q**2
            + specialized_k * (specialized_a2**2 + 2 * d0 * d2)
        )
        * x**2
        + 2 * specialized_k * specialized_a2 * d2 * x**3
        + specialized_k * d2**2 * x**4,
        x,
        domain=sp.QQ,
    )
    specialized_gcd = equation_one.gcd(equation_zero).monic()
    expected_gcd = sp.Poly(x - specialized_b1, x, domain=sp.QQ)
    quotient_resultant = sp.resultant(
        equation_one.exquo(specialized_gcd),
        equation_zero.exquo(specialized_gcd),
        x,
    )
    print("B1_GCD_WITNESS_POINT", (3, 2, 5, 1))
    print("SPECIALIZED_L", sp.Rational(leading_numerator.subs(point)))
    print("SPECIALIZED_K", specialized_k)
    print("SPECIALIZED_B2", specialized_b2)
    print("SPECIALIZED_B1", specialized_b1)
    print("SPECIALIZED_EQUATION_DEGREES", equation_one.degree(), equation_zero.degree())
    print("SPECIALIZED_GCD", specialized_gcd.as_expr())
    print("SPECIALIZED_QUOTIENT_RESULTANT", quotient_resultant)
    if specialized_gcd != expected_gcd or quotient_resultant == 0:
        raise RuntimeError("the exact gcd witness failed")
    print("FUNCTION_FIELD_B1_GCD_DEGREE", 1)

    formulas = {
        name: sp.factor(coefficient_field.to_sympy(value))
        for name, value in (
            ("A2", a2),
            ("A1", a1),
            ("A0", a0),
            ("B2", b2),
            ("B1", b1),
            ("B0", b0),
        )
    }
    lines = []
    for name, expression in formulas.items():
        numerator, denominator = sp.fraction(expression)
        print(
            name,
            "NUM_TERMS",
            len(sp.Poly(numerator, *parameters).terms()),
            "DEN_TERMS",
            len(sp.Poly(denominator, *parameters).terms()),
        )
        lines.append(f"{name}={sp.sstr(expression)}\n")
    formula_text = "".join(lines)
    relative_formula_path = FORMULAS.relative_to(HERE.parent).as_posix()
    if write_formulas:
        FORMULAS.write_text(formula_text, encoding="ascii")
        print("FORMULA_FILE_ACTION", "written")
    else:
        if not FORMULAS.is_file():
            raise RuntimeError(
                f"missing saved formulas: {relative_formula_path}; "
                "rerun with --write-formulas to create them"
            )
        saved_formula_text = FORMULAS.read_text(encoding="ascii")
        if saved_formula_text != formula_text:
            raise RuntimeError(
                "the regenerated formulas differ from the saved certificate; "
                "inspect the inputs before using --write-formulas"
            )
        print("FORMULA_FILE_ACTION", "verified-read-only")
    print("FORMULA_FILE", relative_formula_path)
    print("UNIVERSAL_FIELD_NORM_IDENTITY_REMAINDER", 0)


if __name__ == "__main__":
    arguments = parse_args()
    main(write_formulas=arguments.write_formulas)
