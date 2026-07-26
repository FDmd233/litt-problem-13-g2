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

This script only computes the exact characteristic-zero norm-conic
determinant and emits Singular inputs for factor/genus and ordering-cover
checks.  It does not infer generic behavior from one fiber.
"""

from __future__ import annotations

import pathlib

import sympy as sp


A = -64
B = 64


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
    point = {u: 1, v: 0, z: 9}
    conic_coefficients = [
        sp.expand(coefficient.subs(point))
        for coefficient in data["conic_coefficients"]
    ]
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


def main() -> None:
    data = build_sample()
    print("FIELD QQ")
    print(f"ELLIPTIC_CUBIC {data['elliptic_cubic']}")
    print(f"BRANCH_DISCRIMINANT {data['branch_discriminant']}")
    print(f"FIXED_LINE_INTERSECTION {sp.factor(data['fixed_line_polynomial'])}")
    print(f"FIXED_POINTS {data['fixed_points']}")
    print(f"POINT_CHECKS {data['point_checks']}")
    print(f"EVAL_3K_RANK {data['evaluation_3k'].rank()}")
    print(f"EVAL_2K_RANK {data['evaluation_2k'].rank()}")
    print(f"KERNEL {data['kernel']}")
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
    factor_path, cover_path = write_singular_inputs(data)
    modular_path = write_modular_checks(data)
    print(f"FACTOR_INPUT {factor_path}")
    print(f"COVER_INPUT {cover_path}")
    print(f"MODULAR_INPUT {modular_path}")


if __name__ == "__main__":
    main()
