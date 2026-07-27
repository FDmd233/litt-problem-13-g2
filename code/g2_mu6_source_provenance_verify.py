#!/usr/bin/env python3
"""Regenerate and verify the two universal Singular inputs from the norm.

The calculation starts with the normalized cyclic cover

    E: y^2 = x^2 + lam*(x^3-S1*x^2+S2*x-e3^3),
    C: w^3 = y,

where S1=e1^3-3*e1*e2+3*e3 and
S2=e2^3-3*e1*e2*e3+3*e3^2.  It constructs the three sections vanishing on
the fixed divisor, takes their cyclic cubic norm, divides exactly by y-x,
and computes the residual conic determinant and ordering discriminant.

By default the regenerated sources are compared with the saved files without
writing.  Pass --write-inputs only when intentionally replacing those files.
"""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import re

import sympy as sp


HERE = pathlib.Path(__file__).resolve().parent
REPOSITORY = HERE.parent
GENERATED = HERE / "generated"
E1ONE_VONE_SOURCE = GENERATED / "g2_mu6_universal_cover_e1one_vone.sing"
SPECIAL123_SOURCE = GENERATED / "g2_mu6_universal_cover_special_123.sing"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-inputs",
        action="store_true",
        help="replace the two saved Singular inputs instead of verifying them",
    )
    return parser.parse_args()


def reduce_on_elliptic(
    expression: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    cubic: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    """Return A(x), B(x) with expression=A+y*B modulo y^2=cubic."""
    even = sp.Integer(0)
    odd = sp.Integer(0)
    for (power,), coefficient in sp.Poly(sp.expand(expression), y).terms():
        reduced_coefficient = sp.expand(coefficient * cubic ** (power // 2))
        if power % 2:
            odd += reduced_coefficient
        else:
            even += reduced_coefficient
    return sp.expand(even), sp.expand(odd)


def coefficient_vector(
    expression: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    cubic: sp.Expr,
) -> sp.Matrix:
    """Coordinates in 1,x,...,x^4,y,xy,...,x^3*y in H^0(E,9[O])."""
    even, odd = reduce_on_elliptic(expression, x, y, cubic)
    even_polynomial = sp.Poly(even, x, domain="EX")
    odd_polynomial = sp.Poly(odd, x, domain="EX")
    if even_polynomial.degree() > 4 or odd_polynomial.degree() > 3:
        raise RuntimeError(
            "section outside H^0(E,9[O]): "
            f"degrees {even_polynomial.degree()},{odd_polynomial.degree()}"
        )
    return sp.Matrix(
        [
            *[even_polynomial.nth(degree) for degree in range(5)],
            *[odd_polynomial.nth(degree) for degree in range(4)],
        ]
    )


def singular(expression: sp.Expr) -> str:
    return str(sp.factor(expression)).replace("**", "^")


def read_assignment(
    source: pathlib.Path,
    name: str,
    local_dict: dict[str, sp.Symbol],
) -> sp.Expr:
    text = source.read_text(encoding="ascii")
    match = re.search(rf"^poly {re.escape(name)}=(.*);$", text, re.MULTILINE)
    if match is None:
        relative = source.relative_to(REPOSITORY).as_posix()
        raise RuntimeError(f"missing {name} assignment in {relative}")
    return sp.sympify(match.group(1).replace("^", "**"), locals=local_dict)


def assert_identity(label: str, left: sp.Expr, right: sp.Expr) -> None:
    remainder = sp.cancel(left - right)
    if remainder != 0:
        numerator = sp.factor(sp.fraction(remainder)[0])
        raise RuntimeError(f"{label} failed; numerator remainder: {numerator}")
    print(label, True)


def polynomial_power_mod(
    base: sp.Poly,
    exponent: int,
    modulus: sp.Poly,
) -> sp.Poly:
    """Return base**exponent modulo a monic finite-field polynomial."""
    result = sp.Poly(1, *modulus.gens, modulus=modulus.get_modulus())
    current = base.rem(modulus)
    while exponent:
        if exponent & 1:
            result = (result * current).rem(modulus)
        current = (current * current).rem(modulus)
        exponent >>= 1
    return result


def cover_source(
    *,
    ring: str,
    curve: sp.Expr,
    four_squareclass: sp.Expr,
    list_before_dimension: bool,
) -> str:
    if list_before_dimension:
        decomposition_header = '''list L=minAssGTZ(G);
"ADJUSTED_COVER_VDIM",vdim(G);
"ADJUSTED_COVER_MINASS",size(L);'''
    else:
        decomposition_header = '''"ADJUSTED_COVER_VDIM",vdim(G);
list L=minAssGTZ(G);
"ADJUSTED_COVER_MINASS",size(L);'''
    return f'''LIB "primdec.lib";
ring {ring};
poly curve={singular(curve)};
poly fourS={singular(four_squareclass)};
ideal adjusted=curve,t^2+3*fourS;
ideal G=std(adjusted);
{decomposition_header}
int i;
for(i=1;i<=size(L);i++){{
  ideal H=std(L[i]);
  "ADJUSTED_COVER_ASS",i,vdim(H),size(H);
}}
poly relation=L[1][2];
poly rootnum=subst(relation,t,0);
poly rootden=subst(relation,t,1)-rootnum;
poly identity=reduce(rootnum^2+3*rootden^2*fourS,std(curve));
"ROOT_DENOMINATOR",rootden;
"ROOT_NUMERATOR",rootnum;
"ADJUSTED_SQUARE_IDENTITY_REMAINDER",identity;
quit;
'''


def canonical_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def write_or_verify(
    path: pathlib.Path,
    expected_text: str,
    *,
    write_inputs: bool,
) -> None:
    relative = path.relative_to(REPOSITORY).as_posix()
    if write_inputs:
        path.write_text(expected_text, encoding="ascii")
        action = "written"
    else:
        if not path.is_file():
            raise RuntimeError(
                f"missing generated input {relative}; rerun with --write-inputs"
            )
        actual_text = path.read_text(encoding="ascii")
        if actual_text != expected_text:
            raise RuntimeError(
                f"regenerated source differs from {relative}; "
                "inspect the symbolic identities before using --write-inputs"
            )
        action = "verified-read-only"
    print("SOURCE_FILE_ACTION", relative, action)
    print("SOURCE_CANONICAL_SHA256", relative, canonical_sha256(expected_text))


def main(*, write_inputs: bool = False) -> None:
    x, y, u, v, z = sp.symbols("x y u v z")
    lam, e1, e2, e3 = sp.symbols("lam e1 e2 e3")
    all_symbols = (x, y, u, v, z, lam, e1, e2, e3)
    local_dict = {str(symbol): symbol for symbol in all_symbols}

    sum_cubes = e1**3 - 3 * e1 * e2 + 3 * e3
    pair_cube_sum = e2**3 - 3 * e1 * e2 * e3 + 3 * e3**2
    elliptic_cubic = sp.expand(
        x**2
        + lam
        * (x**3 - sum_cubes * x**2 + pair_cube_sum * x - e3**3)
    )

    # These are the coefficients alpha,beta,gamma of
    # u(x-y)+v(y-e3-e1*w^2+e2*w)+z((x-e3)w-e1*y+e2*w^2).
    alpha = sp.expand(u * (x - y) + v * (y - e3) - e1 * z * y)
    beta = sp.expand(e2 * v + z * (x - e3))
    gamma = sp.expand(-e1 * v + e2 * z)
    norm = sp.expand(
        alpha**3
        + beta**3 * y
        + gamma**3 * y**2
        - 3 * alpha * beta * gamma * y
    )

    conic_basis = [1, x, y, x**2, x * y, y**2]
    multiplication = sp.Matrix.hstack(
        *[
            coefficient_vector(
                (y - x) * monomial,
                x,
                y,
                elliptic_cubic,
            )
            for monomial in conic_basis
        ]
    )
    norm_vector = coefficient_vector(norm, x, y, elliptic_cubic)
    rows = list(range(6))
    square_minor = multiplication[rows, :]
    if sp.factor(square_minor.det()) == 0:
        raise RuntimeError("the fixed norm-division minor is singular")
    quotient = square_minor.inv() * norm_vector[rows, :]
    quotient = sp.Matrix([sp.factor(entry) for entry in quotient])
    residual = multiplication * quotient - norm_vector
    if any(sp.cancel(entry) != 0 for entry in residual):
        raise RuntimeError("exact division of the cyclic norm by y-x failed")
    print("NORM_DIVISION_ROWS", rows)
    print("NORM_DIVISION_REMAINDER", 0)

    q00, q01, q02, q11, q12, q22 = quotient
    conic_matrix = sp.Matrix(
        [
            [q00, q01 / 2, q02 / 2],
            [q01 / 2, q11, q12 / 2],
            [q02 / 2, q12 / 2, q22],
        ]
    )
    determinant = sp.factor(conic_matrix.det())
    delta = sp.factor(-4 * (q00 * q11 - (q01 / 2) ** 2))
    coefficient_field = sp.QQ.frac_field(lam, e1, e2, e3)
    determinant_polynomial = sp.Poly(
        determinant,
        u,
        v,
        z,
        domain=coefficient_field,
    )
    structural_line = sp.Poly(u - v, u, v, z, domain=coefficient_field)
    residual_octic, structural_remainder = determinant_polynomial.div(
        structural_line
    )
    if not structural_remainder.is_zero:
        raise RuntimeError("u-v does not divide the norm-conic determinant")
    if determinant_polynomial.total_degree() != 9 or residual_octic.total_degree() != 8:
        raise RuntimeError("unexpected determinant or residual degree")
    if residual_octic.rem(structural_line).is_zero:
        raise RuntimeError("the structural factor u-v has multiplicity greater than one")
    if any(sum(monomial) != 8 for monomial, _coefficient in residual_octic.terms()):
        raise RuntimeError("the residual factor is not homogeneous of degree eight")
    print("DETERMINANT_DEGREE", determinant_polynomial.total_degree())
    print("STRUCTURAL_FACTOR_REMAINDER", 0)
    print("STRUCTURAL_FACTOR_MULTIPLICITY", 1)
    print("RESIDUAL_FACTOR_DEGREE", residual_octic.total_degree())
    print("RESIDUAL_FACTOR_HOMOGENEOUS", True)
    monic_residual_octic = residual_octic.monic().as_expr()

    e1one_vone_substitution = {e1: 1, v: 1}
    e1one_vone_curve = sp.factor(
        sp.Rational(4, 9)
        * determinant.subs(e1one_vone_substitution)
        / (e2 * (e2 - e3) * (u - 1))
    )
    e1one_vone_delta = sp.factor(delta.subs(e1one_vone_substitution))

    special_substitution = {lam: 1, e1: 6, e2: 11, e3: 6, z: 1}
    special_curve = sp.factor(
        determinant.subs(special_substitution) / (8910 * (u - v))
    )
    special_delta = sp.factor(delta.subs(special_substitution))

    # A small independent irreducibility certificate for the special octic.
    # Clear the constant denominator, reduce modulo 13, set u=1, and normalize.
    special_denominator, special_integer_curve = sp.Poly(
        special_curve, u, v, domain=sp.QQ
    ).clear_denoms(convert=True)
    denominator_mod13 = int(special_denominator) % 13
    if denominator_mod13 == 0:
        raise RuntimeError("the special octic denominator vanishes modulo 13")
    denominator_inverse = pow(denominator_mod13, -1, 13)
    coefficient_ring = sp.GF(13).poly_ring(u)
    special_as_polynomial_in_v = sp.Poly(
        denominator_inverse * special_integer_curve.as_expr(),
        v,
        domain=coefficient_ring,
    )
    if special_as_polynomial_in_v.degree() != 8:
        raise RuntimeError("the special octic does not have full v-degree")
    if special_as_polynomial_in_v.LC() != coefficient_ring.convert(-5):
        raise RuntimeError("unexpected leading v-coefficient in the special octic")

    h = sp.Poly(
        special_integer_curve.as_expr().subs(u, 1), v, modulus=13
    ).monic()
    expected_h = sp.Poly(
        v**8
        + 3 * v**7
        - 4 * v**6
        - 3 * v**5
        - 4 * v**4
        - 5 * v**3
        + 5 * v**2
        - 2 * v
        - 1,
        v,
        modulus=13,
    )
    if h != expected_h:
        raise RuntimeError("the reconstructed Rabin slice differs from the certificate")
    variable = sp.Poly(v, v, modulus=13)
    q4_difference = polynomial_power_mod(variable, 13**4, h) - variable
    q4_gcd = sp.gcd(h, q4_difference).monic()
    q8_remainder = (polynomial_power_mod(variable, 13**8, h) - variable).rem(h)
    if q4_gcd.degree() != 0 or not q8_remainder.is_zero:
        raise RuntimeError("the Rabin irreducibility certificate failed")
    special_mod13 = sp.Poly(
        denominator_inverse * special_integer_curve.as_expr(),
        u,
        v,
        modulus=13,
    )
    witness_substitution = {u: 0, v: 2}
    point_value = int(special_mod13.eval(witness_substitution))
    point_gradient = tuple(
        int(special_mod13.diff(indeterminate).eval(witness_substitution))
        for indeterminate in (u, v)
    )
    if point_value != 0:
        raise RuntimeError("the certified F_13 point does not lie on the special octic")
    if point_gradient == (0, 0):
        raise RuntimeError("the certified F_13 point is singular on the special octic")
    print("RABIN_SLICE_H", h.as_expr())
    print("RABIN_PARENT_V_DEGREE", special_as_polynomial_in_v.degree())
    print("RABIN_PARENT_V_LEADING_COEFFICIENT_MOD13", -5)
    print("RABIN_Q4_GCD", q4_gcd.as_expr())
    print("RABIN_Q8_REMAINDER", 0)
    print("RABIN_SMOOTH_POINT", (0, 2))
    print("RABIN_POINT_VALUE", point_value)
    print("RABIN_POINT_GRADIENT", point_gradient)
    assert_identity(
        "E1ONE_VONE_NORMALIZATION_FROM_MONIC_OCTIC",
        e1one_vone_curve,
        monic_residual_octic.subs(e1one_vone_substitution),
    )
    assert_identity(
        "SPECIAL123_NORMALIZATION_FROM_MONIC_OCTIC",
        special_curve,
        monic_residual_octic.subs(special_substitution),
    )

    def verify_saved_assignments() -> None:
        if not E1ONE_VONE_SOURCE.is_file():
            raise RuntimeError("the saved e1=1,v=1 source is missing")
        saved_e1one_vone_curve = read_assignment(
            E1ONE_VONE_SOURCE, "curve", local_dict
        )
        saved_e1one_vone_delta = read_assignment(
            E1ONE_VONE_SOURCE, "fourS", local_dict
        )
        assert_identity(
            "E1ONE_VONE_CURVE_FROM_DETERMINANT",
            saved_e1one_vone_curve,
            e1one_vone_curve,
        )
        assert_identity(
            "E1ONE_VONE_FOURS_FROM_DELTA",
            saved_e1one_vone_delta,
            e1one_vone_delta,
        )
        if not SPECIAL123_SOURCE.is_file():
            raise RuntimeError("the saved (a,b,c)=(1,2,3) source is missing")
        saved_special_curve = read_assignment(SPECIAL123_SOURCE, "curve", local_dict)
        saved_special_delta = read_assignment(SPECIAL123_SOURCE, "fourS", local_dict)
        assert_identity(
            "SPECIAL123_CURVE_FROM_DETERMINANT",
            saved_special_curve,
            special_curve,
        )
        assert_identity(
            "SPECIAL123_FOURS_FROM_DELTA",
            saved_special_delta,
            special_delta,
        )

    if not write_inputs:
        verify_saved_assignments()

    e1one_vone_text = cover_source(
        ring="r=(0,lam,e2,e3,u),(t,z),lp",
        curve=e1one_vone_curve,
        four_squareclass=e1one_vone_delta,
        list_before_dimension=False,
    )
    special_text = cover_source(
        ring="r=(0,u),(t,v),lp",
        curve=special_curve,
        four_squareclass=special_delta,
        list_before_dimension=True,
    )
    write_or_verify(
        E1ONE_VONE_SOURCE,
        e1one_vone_text,
        write_inputs=write_inputs,
    )
    write_or_verify(
        SPECIAL123_SOURCE,
        special_text,
        write_inputs=write_inputs,
    )
    if write_inputs:
        verify_saved_assignments()


if __name__ == "__main__":
    arguments = parse_args()
    main(write_inputs=arguments.write_inputs)
