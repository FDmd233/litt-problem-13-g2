$ErrorActionPreference = "Stop"

$generated = Join-Path $PSScriptRoot "generated"
$source = Join-Path $generated "g2_mu6_universal_cover_special_123.sing"
$text = Get-Content $source -Raw
$match = [regex]::Match(
    $text,
    'poly curve=\((.*?)\)/35640;',
    [System.Text.RegularExpressions.RegexOptions]::Singleline
)
if (-not $match.Success) {
    throw "Could not extract the special (a,b,c)=(1,2,3) octic."
}
$curve = $match.Groups[1].Value

$singular = @"
ring r=13,(u,v),dp;
poly curve=$curve;
list F=factorize(curve);
"FACTOR_LIST_SIZE",size(F[1]);
"UNIT_FACTOR_DEGREE",deg(F[1][1]);
"NONUNIT_FACTOR_COUNT",size(F[1])-1;
"NONUNIT_FACTOR_DEGREE",deg(F[1][2]);
"NONUNIT_FACTOR_EXPONENT",F[2][2];
poly du=diff(curve,u);
poly dv=diff(curve,v);
"POINT_VALUE",subst(subst(curve,u,0),v,2);
"POINT_GRADIENT",
  subst(subst(du,u,0),v,2),
  subst(subst(dv,u,0),v,2);
quit;
"@

$output = @($singular | wsl /usr/bin/Singular -q 2>&1)
$singularExitCode = $LASTEXITCODE
$output | Write-Output

if ($singularExitCode -ne 0) {
    throw "Singular exited with code $singularExitCode."
}

$joinedOutput = [string]::Join(
    "`n",
    @($output | ForEach-Object { [string]$_ })
)

function Assert-IntegerField {
    param(
        [Parameter(Mandatory = $true)][string]$Name,
        [Parameter(Mandatory = $true)][int]$Expected
    )

    $pattern = "(?m)^" + [regex]::Escape($Name) + "\s+(-?\d+)\s*$"
    $fieldMatches = [regex]::Matches($joinedOutput, $pattern)
    if ($fieldMatches.Count -ne 1) {
        throw "Expected exactly one integer certificate field $Name; found $($fieldMatches.Count)."
    }
    $fieldMatch = $fieldMatches[0]
    $actual = [int]$fieldMatch.Groups[1].Value
    if ($actual -ne $Expected) {
        throw "Certificate field $Name is $actual; expected $Expected."
    }
}

Assert-IntegerField "FACTOR_LIST_SIZE" 2
Assert-IntegerField "UNIT_FACTOR_DEGREE" 0
Assert-IntegerField "NONUNIT_FACTOR_COUNT" 1
Assert-IntegerField "NONUNIT_FACTOR_DEGREE" 8
Assert-IntegerField "NONUNIT_FACTOR_EXPONENT" 1
Assert-IntegerField "POINT_VALUE" 0

$gradientMatches = [regex]::Matches(
    $joinedOutput,
    '(?m)^POINT_GRADIENT\s+(-?\d+)\s+(-?\d+)\s*$'
)
if ($gradientMatches.Count -ne 1) {
    throw "Expected exactly one two-coordinate POINT_GRADIENT field; found $($gradientMatches.Count)."
}
$gradientMatch = $gradientMatches[0]
$gradientU = [int]$gradientMatch.Groups[1].Value
$gradientV = [int]$gradientMatch.Groups[2].Value
$gradientUMod13 = (($gradientU % 13) + 13) % 13
$gradientVMod13 = (($gradientV % 13) + 13) % 13
if (($gradientUMod13 -eq 0) -and ($gradientVMod13 -eq 0)) {
    throw "The certified point is singular: POINT_GRADIENT vanishes modulo 13."
}

Write-Output "CERTIFICATE_ASSERTIONS_OK True"
