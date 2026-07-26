#!/usr/bin/env python3
"""Independently verify the cleared-numerator norm identity.

This script reads the six saved coefficient formulas and the actively generated
curve source.  It does not re-run the recovery equations.  With L the leading
numerator of the source sextic, it forms polynomial numerators

    A_num = 2*L*A,  B_num = 2*L*B,
    C_num = L*monic(curve),  S_num = lam^2*fourS,

and checks in Q[z,lam,e2,e3,u] that

    lam^2*A_num^2 + 3*S_num*B_num^2 - 4*lam^2*L*C_num = 0.
"""

from __future__ import annotations

import pathlib
import re

import sympy as sp


HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "generated" / "g2_mu6_universal_cover_e1one_vone.sing"
FORMULAS = HERE / "generated" / "g2_mu6_e1one_vone_norm_formulas.txt"


def read_source_assignment(name: str) -> str:
    text = SOURCE.read_text(encoding="ascii")
    match = re.search(rf"^poly {re.escape(name)}=(.*);$", text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"missing {name} assignment in {SOURCE}")
    return match.group(1).replace("^", "**")


def polynomial_numerator(label: str, expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    if denominator != 1:
        raise RuntimeError(f"{label} retained denominator {denominator}")
    return sp.expand(numerator)


def main() -> None:
    lam, e2, e3, u, z = sp.symbols("lam e2 e3 u z")
    local_dict = {str(symbol): symbol for symbol in (lam, e2, e3, u, z)}
    curve = sp.sympify(read_source_assignment("curve"), locals=local_dict)
    four_s = sp.sympify(read_source_assignment("fourS"), locals=local_dict)

    formulas: dict[str, sp.Expr] = {}
    for line in FORMULAS.read_text(encoding="ascii").splitlines():
        name, value = line.split("=", 1)
        formulas[name] = sp.sympify(value, locals=local_dict)
    expected_names = {"A2", "A1", "A0", "B2", "B1", "B0"}
    if set(formulas) != expected_names:
        raise RuntimeError(f"bad formula names: {sorted(formulas)}")

    curve_polynomial = sp.Poly(curve, z, domain="EX")
    leading = sp.cancel(curve_polynomial.LC())
    leading_numerator = sp.fraction(leading)[0]
    monic_curve = sp.cancel(curve / leading)
    a = z**3 + formulas["A2"] * z**2 + formulas["A1"] * z + formulas["A0"]
    b = formulas["B2"] * z**2 + formulas["B1"] * z + formulas["B0"]

    a_numerator = polynomial_numerator("2*L*A", 2 * leading_numerator * a)
    b_numerator = polynomial_numerator("2*L*B", 2 * leading_numerator * b)
    curve_numerator = polynomial_numerator("L*monic(curve)", leading_numerator * monic_curve)
    four_s_numerator = polynomial_numerator("lam^2*fourS", lam**2 * four_s)

    generators = (z, lam, e2, e3, u)
    lam_polynomial = sp.Poly(lam, *generators, domain=sp.QQ)
    leading_polynomial = sp.Poly(leading_numerator, *generators, domain=sp.QQ)
    a_polynomial = sp.Poly(a_numerator, *generators, domain=sp.QQ)
    b_polynomial = sp.Poly(b_numerator, *generators, domain=sp.QQ)
    curve_numerator_polynomial = sp.Poly(curve_numerator, *generators, domain=sp.QQ)
    four_s_numerator_polynomial = sp.Poly(four_s_numerator, *generators, domain=sp.QQ)

    print("A_NUMERATOR_TERMS", len(a_polynomial.terms()))
    print("B_NUMERATOR_TERMS", len(b_polynomial.terms()))
    print("CURVE_NUMERATOR_TERMS", len(curve_numerator_polynomial.terms()))
    print("FOURS_NUMERATOR_TERMS", len(four_s_numerator_polynomial.terms()))
    certificate = (
        lam_polynomial**2 * a_polynomial**2
        + 3 * four_s_numerator_polynomial * b_polynomial**2
        - 4
        * lam_polynomial**2
        * leading_polynomial
        * curve_numerator_polynomial
    )
    print("CLEARED_COMMON_DENOMINATOR", "4*lam^2*L^2")
    print("DIRECT_NUMERATOR_TERMS", 0 if certificate.is_zero else len(certificate.terms()))
    print("DIRECT_NUMERATOR_IS_ZERO", certificate.is_zero)
    if not certificate.is_zero:
        raise RuntimeError("the direct cleared-numerator identity failed")
    print("DIRECT_UNIVERSAL_NORM_IDENTITY_REMAINDER", 0)


if __name__ == "__main__":
    main()
