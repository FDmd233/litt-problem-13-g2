#!/usr/bin/env python3
"""Exact norm-conic computation for a nonsimple mu_6-symmetric G2 sample.

The sample is

    E: y^2 = x^3 - 64*x + 64,
    C: w^3 = y,

so equivalently C is the superelliptic curve w^6=x^3-64*x+64.  Besides
the cyclic deck transformation w -> zeta_3*w, it has the involution
w -> -w.  The latter has elliptic quotient u^3=x^3-64*x+64 of j-invariant
zero, and the Prym is nonsimple.

The fixed divisor is cut out by y=x and consists of

    (1,1,1), (8,8,2), (-8,-8,-2).

This script computes the exact characteristic-zero norm-conic determinant.
It also recomputes, with exact SymPy arithmetic, every modulo-13 value in the
finite-field open witness printed in the paper.  The emitted Singular inputs
are optional cross-checks; they are not needed for the Python assertions.
No generic behavior is inferred from one fiber.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable
import pathlib

import sympy as sp


A = -64
B = 64
OPEN_WITNESS_PRIME = 13


def reduce_on_elliptic(
    expression: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> tuple[sp.Expr, sp.Expr]:
    """Return A(x),B(x) with expression=A+y*B modulo the elliptic equation."""
    cubic = x**3 + A * x + B
    even = sp.Integer(0)
    odd = sp.Integer(0)
    for (power,), coefficient in sp.Poly(sp.expand(expression), y).terms():
        target = sp.expand(coefficient * cubic ** (power // 2))
        if power % 2:
            odd += target
        else:
            even += target
    return sp.expand(even), sp.expand(odd)


def coefficient_vector(
    expression: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> sp.Matrix:
    """Coordinates in 1,x,...,x^4,y,xy,...,x^3*y."""
    even, odd = reduce_on_elliptic(expression, x, y)
    even_poly = sp.Poly(even, x, domain="EX")
    odd_poly = sp.Poly(odd, x, domain="EX")
    if even_poly.degree() > 4 or odd_poly.degree() > 3:
        raise RuntimeError(
            f"section is outside H^0(E,9[O]): degrees "
            f"{even_poly.degree()},{odd_poly.degree()}"
        )
    return sp.Matrix(
        [
            *[even_poly.nth(i) for i in range(5)],
            *[odd_poly.nth(i) for i in range(4)],
        ]
    )


def independent_rows(matrix: sp.Matrix) -> list[int]:
    chosen: list[int] = []
    rank = 0
    for row in range(matrix.rows):
        trial_rank = matrix[chosen + [row], :].rank()
        if trial_rank > rank:
            chosen.append(row)
            rank = trial_rank
    return chosen


def primitive_integer_vector(vector: sp.Matrix) -> list[int]:
    denominators = [sp.denom(entry) for entry in vector]
    common_denominator = sp.ilcm(*[int(value) for value in denominators])
    values = [int(entry * common_denominator) for entry in vector]
    common_divisor = abs(sp.igcd(*values))
    values = [value // common_divisor for value in values]
    if next(value for value in values if value) < 0:
        values = [-value for value in values]
    return values


def singular(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def residue_mod_prime(value: sp.Expr, prime: int) -> int:
    """Reduce an exact rational value modulo prime."""
    rational = sp.Rational(sp.cancel(value))
    numerator = int(sp.numer(rational))
    denominator = int(sp.denom(rational))
    if denominator % prime == 0:
        raise ZeroDivisionError(
            f"denominator {denominator} is zero modulo {prime}"
        )
    return numerator * pow(denominator, -1, prime) % prime


def residue_tuple(values: Iterable[sp.Expr], prime: int) -> tuple[int, ...]:
    return tuple(residue_mod_prime(value, prime) for value in values)


def assert_equal(label: str, actual: object, expected: object) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected}, got {actual}")


def build_sample() -> dict[str, object]:
    x, y, u, v, z = sp.symbols("x y u v z")
    parameters = (u, v, z)
    fixed_points = [(1, 1, 1), (8, 8, 2), (-8, -8, -2)]

    elliptic_cubic = x**3 + A * x + B
    branch_polynomial = elliptic_cubic
    fixed_line_polynomial = sp.expand(x**2 - elliptic_cubic)
    point_checks = [
        (
            y_value**2 - (x_value**3 + A * x_value + B),
            w_value**3 - y_value,
            y_value - x_value,
            y_value,
        )
        for x_value, y_value, w_value in fixed_points
    ]

    evaluation_3k = sp.Matrix(
        [
            [1, x_value, y_value, w_value, x_value * w_value, w_value**2]
            for x_value, y_value, w_value in fixed_points
        ]
    )
    evaluation_2k = sp.Matrix(
        [
            [1, x_value, w_value, w_value**2]
            for x_value, _y_value, w_value in fixed_points
        ]
    )
    kernel = [primitive_integer_vector(vector) for vector in evaluation_3k.nullspace()]
    if evaluation_3k.rank() != 3 or len(kernel) != 3:
        raise RuntimeError("evaluation on H^0(3*kappa) does not have rank three")
    if evaluation_2k.rank() != 3:
        raise RuntimeError("fixed divisor does not map to a smooth theta point")

    coefficients = [
        sp.expand(sum(parameters[i] * kernel[i][j] for i in range(3)))
        for j in range(6)
    ]
    aa = coefficients[0] + coefficients[1] * x + coefficients[2] * y
    bb = coefficients[3] + coefficients[4] * x
    cc = coefficients[5]
    norm = sp.expand(aa**3 + bb**3 * y + cc**3 * y**2 - 3 * aa * bb * cc * y)

    fixed_line = y - x
    conic_basis = [1, x, y, x**2, x * y, y**2]
    multiplication_matrix = sp.Matrix.hstack(
        *[coefficient_vector(fixed_line * monomial, x, y) for monomial in conic_basis]
    )
    norm_vector = coefficient_vector(norm, x, y)
    rows = independent_rows(multiplication_matrix)
    if len(rows) != 6:
        raise RuntimeError("multiplication by the fixed line is not injective")
    quotient = multiplication_matrix[rows, :].inv() * norm_vector[rows, :]
    if sp.simplify(multiplication_matrix * quotient - norm_vector) != sp.zeros(9, 1):
        raise RuntimeError("norm division failed")

    q00, q01, q02, q11, q12, q22 = map(sp.expand, quotient)
    conic_matrix = sp.Matrix(
        [
            [q00, q01 / 2, q02 / 2],
            [q01 / 2, q11, q12 / 2],
            [q02 / 2, q12 / 2, q22],
        ]
    )
    determinant = sp.cancel(sp.expand(conic_matrix.det()))
    determinant_numerator, determinant_denominator = sp.fraction(determinant)
    determinant_poly = sp.Poly(
        determinant_numerator, u, v, z, domain=sp.QQ
    ).primitive()[1]
    factorization = sp.factor_list(determinant_poly.as_expr(), u, v, z)
    squareclass = sp.Poly(
        -(q00 * q11 - (q01 / 2) ** 2), u, v, z, domain=sp.QQ
    )

    return {
        "symbols": (x, y, u, v, z),
        "elliptic_cubic": elliptic_cubic,
        "branch_discriminant": sp.discriminant(branch_polynomial, x),
        "fixed_line_polynomial": fixed_line_polynomial,
        "fixed_points": fixed_points,
        "point_checks": point_checks,
        "evaluation_3k": evaluation_3k,
        "evaluation_2k": evaluation_2k,
        "kernel": kernel,
        "conic_coefficients": tuple(map(sp.expand, quotient)),
        "multiplication_rank": multiplication_matrix.rank(),
        "determinant_denominator": determinant_denominator,
        "determinant": determinant_poly,
        "factorization": factorization,
        "squareclass": squareclass,
    }


def open_witness(data: dict[str, object]) -> dict[str, object]:
    """Express the finite-field open witness in both coordinate bases."""
    _x, _y, u, v, z = data["symbols"]
    parameters = (u, v, z)
    kernel_point = {u: 1, v: 0, z: 9}
    kernel_section = sp.Matrix(
        [
            sp.expand(
                sum(parameters[i] * data["kernel"][i][j] for i in range(3))
            ).subs(kernel_point)
            for j in range(6)
        ]
    )

    # Coefficient order: 1,x,y,w,xw,w^2.  These are the three sections
    # s_0,s_1,s_2 displayed in the paper for (a,b,c)=(1,2,-2).
    paper_basis = sp.Matrix.hstack(
        sp.Matrix([0, 1, -1, 0, 0, 0]),
        sp.Matrix([4, 0, 1, -4, 0, -1]),
        sp.Matrix([0, 0, -1, 4, 1, -4]),
    )
    coordinate_solutions = sp.linsolve((paper_basis, kernel_section))
    if len(coordinate_solutions) != 1:
        raise RuntimeError("the witness section has no unique paper-basis coordinates")
    paper_point = tuple(next(iter(coordinate_solutions)))
    if paper_point != (sp.Integer(10), sp.Integer(9), sp.Integer(0)):
        raise RuntimeError(f"unexpected paper-basis witness point: {paper_point}")

    kernel_to_paper = sp.Matrix(
        [
            [1, 5, 1],
            [0, 4, 1],
            [0, -1, 0],
        ]
    )
    if paper_basis * kernel_to_paper != sp.Matrix.hstack(
        *[sp.Matrix(vector) for vector in data["kernel"]]
    ):
        raise RuntimeError("the kernel-to-paper basis conversion failed")

    conic_coefficients = [
        sp.expand(coefficient.subs(kernel_point))
        for coefficient in data["conic_coefficients"]
    ]
    q00, q01, q02, q11, q12, q22 = conic_coefficients
    conic_matrix = sp.Matrix(
        [
            [q00, q01 / 2, q02 / 2],
            [q01 / 2, q11, q12 / 2],
            [q02 / 2, q12 / 2, q22],
        ]
    )
    conic_determinant = sp.factor(conic_matrix.det())
    rank_two_minor = sp.expand(4 * q00 * q11 - q01**2)
    determinant_mod13 = residue_mod_prime(conic_determinant, OPEN_WITNESS_PRIME)
    rank_two_minor_mod13 = residue_mod_prime(rank_two_minor, OPEN_WITNESS_PRIME)
    assert_equal("universal-basis witness point", paper_point, (10, 9, 0))
    assert_equal("witness conic determinant modulo 13", determinant_mod13, 0)
    assert_equal("witness rank-two minor modulo 13", rank_two_minor_mod13, 10)

    return {
        "kernel_point": (1, 0, 9),
        "paper_point": tuple(map(int, paper_point)),
        "kernel_to_paper": kernel_to_paper,
        "section_vector": tuple(map(int, kernel_section)),
        "conic_coefficients": conic_coefficients,
        "conic_determinant": conic_determinant,
        "rank_two_minor": rank_two_minor,
        "determinant_mod13": determinant_mod13,
        "rank_two_minor_mod13": rank_two_minor_mod13,
    }


def finite_field_open_checks(
    data: dict[str, object], witness: dict[str, object]
) -> dict[str, object]:
    """Recompute and assert the complete F_13 open-locus witness."""
    x, y, _u, _v, _z = data["symbols"]
    prime = OPEN_WITNESS_PRIME

    branch_discriminant_mod13 = residue_mod_prime(
        data["branch_discriminant"], prime
    )
    evaluation_3k_minor_mod13 = residue_mod_prime(
        data["evaluation_3k"].extract([0, 1, 2], [0, 1, 3]).det(), prime
    )
    evaluation_2k_minor_mod13 = residue_mod_prime(
        data["evaluation_2k"].extract([0, 1, 2], [0, 1, 2]).det(), prime
    )
    assert_equal(
        "elliptic branch discriminant modulo 13", branch_discriminant_mod13, 8
    )
    assert_equal(
        "3-kappa evaluation minor modulo 13", evaluation_3k_minor_mod13, 1
    )
    assert_equal(
        "2-kappa evaluation minor modulo 13", evaluation_2k_minor_mod13, 1
    )

    conic_coefficients_mod13 = residue_tuple(
        witness["conic_coefficients"], prime
    )
    assert_equal(
        "witness conic coefficients modulo 13",
        conic_coefficients_mod13,
        (4, 3, 1, 2, 3, 12),
    )

    q00, q01, q02, q11, q12, q22 = witness["conic_coefficients"]
    conic = sp.expand(
        q00 + q01 * x + q02 * y + q11 * x**2 + q12 * x * y + q22 * y**2
    )
    line1 = 4 * x + y - 6
    line2 = 6 * x + y + 5
    factor_remainder = sp.Poly(conic + line1 * line2, x, y, modulus=prime)
    assert_equal(
        "witness conic factorization modulo 13", factor_remainder.is_zero, True
    )

    fixed_point_line_values_exact = tuple(
        line.subs({x: x_value, y: y_value})
        for line in (line1, line2)
        for x_value, y_value, _w_value in data["fixed_points"]
    )
    fixed_point_line_values = residue_tuple(fixed_point_line_values_exact, prime)
    assert_equal(
        "fixed-point line values modulo 13",
        fixed_point_line_values,
        (12, 8, 6, 12, 9, 1),
    )

    elliptic_cubic = data["elliptic_cubic"]
    line_y_coordinates = (6 - 4 * x, -6 * x - 5)
    line_sections = tuple(
        sp.expand(y_coordinate**2 - elliptic_cubic)
        for y_coordinate in line_y_coordinates
    )
    derivative_resultants_exact = tuple(
        sp.resultant(section, sp.diff(section, x), x)
        for section in line_sections
    )
    branch_resultants_exact = tuple(
        sp.resultant(section, y_coordinate, x)
        for section, y_coordinate in zip(line_sections, line_y_coordinates)
    )
    pair_resultant_exact = sp.resultant(line_sections[0], line_sections[1], x)
    derivative_resultants = residue_tuple(derivative_resultants_exact, prime)
    branch_resultants = residue_tuple(branch_resultants_exact, prime)
    pair_resultant = residue_mod_prime(pair_resultant_exact, prime)
    assert_equal(
        "line-section derivative resultants modulo 13",
        derivative_resultants,
        (10, 6),
    )
    assert_equal(
        "line-branch resultants modulo 13", branch_resultants, (12, 1)
    )
    assert_equal("line-pair resultant modulo 13", pair_resultant, 12)

    return {
        "prime": prime,
        "branch_discriminant_mod13": branch_discriminant_mod13,
        "evaluation_3k_minor_mod13": evaluation_3k_minor_mod13,
        "evaluation_2k_minor_mod13": evaluation_2k_minor_mod13,
        "conic_coefficients_mod13": conic_coefficients_mod13,
        "conic_factor_remainder_mod13": factor_remainder.as_expr(),
        "fixed_point_line_values_exact": fixed_point_line_values_exact,
        "fixed_point_line_values": fixed_point_line_values,
        "derivative_resultants_exact": derivative_resultants_exact,
        "derivative_resultants": derivative_resultants,
        "branch_resultants_exact": branch_resultants_exact,
        "branch_resultants": branch_resultants,
        "pair_resultant_exact": pair_resultant_exact,
        "pair_resultant": pair_resultant,
    }


def write_singular_inputs(data: dict[str, object]) -> tuple[pathlib.Path, pathlib.Path]:
    _x, _y, u, v, z = data["symbols"]
    generated = pathlib.Path(__file__).resolve().parent / "generated"
    generated.mkdir(parents=True, exist_ok=True)
    determinant = data["determinant"].as_expr()
    factors = []
    for expression, exponent in data["factorization"][1]:
        factor = sp.Poly(expression, u, v, z, domain=sp.QQ).primitive()[1]
        factors.extend([factor] * exponent)

    factor_path = generated / "g2_mu6_nonsimple_norm_conic_Q_factor.sing"
    declarations = "\n".join(
        f"poly f{index}={singular(factor.as_expr())};"
        for index, factor in enumerate(factors, start=1)
    )
    audits = "\n".join(
        f'"FACTOR",{index},"DEGREE",deg(f{index}),"GENUS",genus(ideal(f{index}));\n'
        f"factorize(f{index});"
        for index in range(1, len(factors) + 1)
    )
    content = data["factorization"][0]
    product = "*".join(f"f{index}" for index in range(1, len(factors) + 1))
    factor_path.write_text(
        f'''LIB "normal.lib";
ring r=0,(u,v,z),dp;
poly determinant={singular(determinant)};
{declarations}
"DETERMINANT_DEGREE",deg(determinant);
"PRODUCT_REMAINDER",determinant-({singular(content)})*({product});
{audits}
quit;
''',
        encoding="ascii",
    )

    cover_path = generated / "g2_mu6_nonsimple_norm_conic_Q_covers.sing"
    cover_blocks = []
    four_squareclass = 4 * data["squareclass"].as_expr().subs(z, 1)
    for index, factor in enumerate(factors, start=1):
        if factor.total_degree() <= 1:
            continue
        affine_factor = factor.as_expr().subs(z, 1)
        cover_blocks.append(
            f'''poly curve{index}={singular(affine_factor)};
poly fourS{index}={singular(four_squareclass)};
ideal base{index}=std(curve{index});
ideal cover{index}=curve{index},t^2-fourS{index};
ideal G{index}=std(cover{index});
list L{index}=minAssGTZ(G{index});
"ARITHMETIC_COVER",{index},"VDIM",vdim(G{index}),"MINASS",size(L{index});
int i{index};
for(i{index}=1;i{index}<=size(L{index});i{index}++){{
  ideal H=std(L{index}[i{index}]);
  "ARITHMETIC_COVER_ASS",{index},i{index},vdim(H),size(H);
}}
ideal adjusted{index}=curve{index},t^2+3*fourS{index};
ideal adjustedG{index}=std(adjusted{index});
list adjustedL{index}=minAssGTZ(adjustedG{index});
"ADJUSTED_COVER",{index},"VDIM",vdim(adjustedG{index}),
  "MINASS",size(adjustedL{index});
for(i{index}=1;i{index}<=size(adjustedL{index});i{index}++){{
  ideal H=std(adjustedL{index}[i{index}]);
  "ADJUSTED_COVER_ASS",{index},i{index},vdim(H),size(H);
}}
poly relation{index}=adjustedL{index}[1][2];
poly rootnum{index}=subst(relation{index},t,0);
poly rootden{index}=subst(relation{index},t,1)-rootnum{index};
poly identity{index}=reduce(
  rootnum{index}^2+3*rootden{index}^2*fourS{index},base{index});
"ROOT_DENOMINATOR",{index},rootden{index};
"ROOT_NUMERATOR",{index},rootnum{index};
"ADJUSTED_SQUARE_IDENTITY_REMAINDER",{index},identity{index};
'''
        )
    cover_path.write_text(
        '''LIB "primdec.lib";
ring r=(0,u),(t,v),lp;
''' + "\n".join(cover_blocks) + "\nquit;\n",
        encoding="ascii",
    )
    return factor_path, cover_path


def write_modular_checks(data: dict[str, object]) -> pathlib.Path:
    """Emit absolute-irreducibility and genuine-open finite-field witnesses."""
    x, y, u, v, z = data["symbols"]
    octics = [
        sp.Poly(expression, u, v, z, domain=sp.QQ).primitive()[1]
        for expression, exponent in data["factorization"][1]
        for _unused in range(exponent)
        if sp.Poly(expression, u, v, z, domain=sp.QQ).total_degree() == 8
    ]
    if len(octics) != 1:
        raise RuntimeError("expected one residual octic")
    octic = octics[0].as_expr()
    witness = open_witness(data)
    conic_coefficients = witness["conic_coefficients"]
    q00, q01, q02, q11, q12, q22 = conic_coefficients
    conic = q00 + q01 * x + q02 * y + q11 * x**2 + q12 * x * y + q22 * y**2

    generated = pathlib.Path(__file__).resolve().parent / "generated"
    modular_path = generated / "g2_mu6_nonsimple_norm_conic_Q_modular.sing"
    modular_path.write_text(
        f'''ring geometric=5,(u,v,z),dp;
poly octic={singular(octic)};
poly du=diff(octic,u);
poly dv=diff(octic,v);
poly dz=diff(octic,z);
"OCTIC_FACTORIZATION_MOD5";
factorize(octic);
"SMOOTH_POINT_VALUE",subst(subst(subst(octic,u,1),v,0),z,2);
"SMOOTH_POINT_GRADIENT",
  subst(subst(subst(du,u,1),v,0),z,2),
  subst(subst(subst(dv,u,1),v,0),z,2),
  subst(subst(subst(dz,u,1),v,0),z,2);

ring opencheck=13,(x,y),dp;
poly conic={singular(conic)};
poly line1=4*x+y-6;
poly line2=6*x+y+5;
"KERNEL_BASIS_PARAMETER_POINT",{','.join(map(str, witness['kernel_point']))};
"PAPER_BASIS_PARAMETER_POINT",{','.join(map(str, witness['paper_point']))};
"CONIC_DETERMINANT",{witness['determinant_mod13']};
"RANK_TWO_MINOR",{witness['rank_two_minor_mod13']};
"PARAMETER_POINT_CONIC_COEFFICIENTS",
  {','.join(singular(entry) for entry in conic_coefficients)};
"CONIC_FACTOR_IDENTITY",conic+line1*line2;
"FIXED_POINT_LINE_VALUES",
  4*1+1-6,4*8+8-6,4*5+5-6,
  6*1+1+5,6*8+8+5,6*5+5+5;

ring resultants=13,(x),dp;
poly elliptic=x^3-64*x+64;
poly y1=6-4*x;
poly y2=-6*x-5;
poly section1=y1^2-elliptic;
poly section2=y2^2-elliptic;
"LINE_SECTION_DERIVATIVE_RESULTANTS",
  resultant(section1,diff(section1,x),x),
  resultant(section2,diff(section2,x),x);
"LINE_BRANCH_RESULTANTS",
  resultant(section1,y1,x),resultant(section2,y2,x);
"LINE_PAIR_RESULTANT",resultant(section1,section2,x);
quit;
''',
        encoding="ascii",
    )
    return modular_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify the exact norm-conic and finite-field open witness."
    )
    parser.add_argument(
        "--write-inputs",
        action="store_true",
        help="regenerate the optional Singular input files",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = build_sample()
    witness = open_witness(data)
    open_checks = finite_field_open_checks(data, witness)
    print("FIELD QQ")
    print(f"ELLIPTIC_CUBIC {data['elliptic_cubic']}")
    print(f"BRANCH_DISCRIMINANT {data['branch_discriminant']}")
    print(f"FIXED_LINE_INTERSECTION {sp.factor(data['fixed_line_polynomial'])}")
    print(f"FIXED_POINTS {data['fixed_points']}")
    print(f"POINT_CHECKS {data['point_checks']}")
    print(f"EVAL_3K_RANK {data['evaluation_3k'].rank()}")
    print(f"EVAL_2K_RANK {data['evaluation_2k'].rank()}")
    print(f"KERNEL {data['kernel']}")
    print(f"KERNEL_TO_PAPER_MATRIX {witness['kernel_to_paper'].tolist()}")
    print(f"WITNESS_KERNEL_BASIS_POINT {witness['kernel_point']}")
    print(f"WITNESS_PAPER_BASIS_POINT {witness['paper_point']}")
    print(f"WITNESS_SECTION_VECTOR {witness['section_vector']}")
    print(
        "BRANCH_DISCRIMINANT_MOD13 "
        f"{open_checks['branch_discriminant_mod13']}"
    )
    print(
        "EVAL_3K_MINOR_COLUMNS_0_1_3_MOD13 "
        f"{open_checks['evaluation_3k_minor_mod13']}"
    )
    print(
        "EVAL_2K_MINOR_COLUMNS_0_1_2_MOD13 "
        f"{open_checks['evaluation_2k_minor_mod13']}"
    )
    print(
        "WITNESS_CONIC_COEFFICIENTS_MOD13 "
        f"{open_checks['conic_coefficients_mod13']}"
    )
    print(f"WITNESS_CONIC_DETERMINANT_MOD13 {witness['determinant_mod13']}")
    print(f"WITNESS_RANK_TWO_MINOR_MOD13 {witness['rank_two_minor_mod13']}")
    print(
        "WITNESS_CONIC_FACTOR_REMAINDER_MOD13 "
        f"{open_checks['conic_factor_remainder_mod13']}"
    )
    print(f"FIXED_POINT_LINE_VALUES {open_checks['fixed_point_line_values']}")
    print(
        "LINE_SECTION_DERIVATIVE_RESULTANTS "
        f"{open_checks['derivative_resultants']}"
    )
    print(f"LINE_BRANCH_RESULTANTS {open_checks['branch_resultants']}")
    print(f"LINE_PAIR_RESULTANT {open_checks['pair_resultant']}")
    print(f"MULTIPLICATION_BY_LINE_RANK {data['multiplication_rank']}")
    print(f"DETERMINANT_DENOMINATOR {data['determinant_denominator']}")
    print(f"DETERMINANT_DEGREE {data['determinant'].total_degree()}")
    print(f"DETERMINANT_TERMS {len(data['determinant'].terms())}")
    print("FACTORIZATION")
    for expression, exponent in data["factorization"][1]:
        polynomial = sp.Poly(expression, *data["symbols"][2:], domain=sp.QQ)
        print(
            f"DEGREE={polynomial.total_degree()} EXPONENT={exponent} "
            f"TERMS={len(polynomial.terms())} FACTOR={polynomial.as_expr()}"
        )
    if args.write_inputs:
        factor_path, cover_path = write_singular_inputs(data)
        modular_path = write_modular_checks(data)
        print(f"FACTOR_INPUT {factor_path}")
        print(f"COVER_INPUT {cover_path}")
        print(f"MODULAR_INPUT {modular_path}")


if __name__ == "__main__":
    main()
