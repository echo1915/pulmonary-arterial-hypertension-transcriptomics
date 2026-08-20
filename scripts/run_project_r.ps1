param(
    [string]$Expression,
    [string]$Script,
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$RArguments
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$RExecutable = "C:\R\R-4.6.1\bin\Rscript.exe"

if (-not (Test-Path -LiteralPath $RExecutable)) {
    throw "R was not found at $RExecutable"
}

$env:R_LIBS_USER = "G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics\cache\R\4.6-library"
$env:TMPDIR = Join-Path $ProjectRoot ".tmp\r-temp"
$env:TEMP = $env:TMPDIR
$env:TMP = $env:TMPDIR
$env:LC_ALL = "C"
$env:LANG = "C"

New-Item -ItemType Directory -Force -Path $env:TMPDIR | Out-Null
if ($Expression) {
    & $RExecutable --vanilla -e $Expression @RArguments
} elseif ($Script) {
    & $RExecutable --vanilla $Script @RArguments
} else {
    & $RExecutable --vanilla @RArguments
}
exit $LASTEXITCODE
