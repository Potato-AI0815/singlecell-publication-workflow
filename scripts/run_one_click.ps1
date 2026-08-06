[CmdletBinding()]
param(
  [Parameter(Mandatory = $true)]
  [string]$Config,

  [ValidateSet('full','fast','resume','reannotate','figures_only','audit','project_profile')]
  [string]$Mode = 'full',

  [string]$RscriptPath,
  [string]$RShortcutPath,
  [switch]$InstallPackages,
  [switch]$InstallGitHubPackages,
  [switch]$RepairCopykat,
  [switch]$SkipDependencyCheck
)

$ErrorActionPreference = 'Stop'
$SkillRoot = Split-Path -Parent $PSScriptRoot
$Config = [System.IO.Path]::GetFullPath($Config)
if (-not (Test-Path -LiteralPath $Config -PathType Leaf)) {
  throw "Configuration file not found: $Config"
}

function Resolve-Rscript {
  param([string]$ExplicitRscript, [string]$ExplicitRShortcut)
  if ($ExplicitRscript) {
    $p = [System.IO.Path]::GetFullPath($ExplicitRscript)
    if (-not (Test-Path -LiteralPath $p -PathType Leaf)) { throw "Rscript not found: $p" }
    return $p
  }
  if ($ExplicitRShortcut) {
    $shortcut = [System.IO.Path]::GetFullPath($ExplicitRShortcut)
    if (-not (Test-Path -LiteralPath $shortcut -PathType Leaf)) { throw "R executable/shortcut not found: $shortcut" }
    $candidate = Join-Path (Split-Path -Parent $shortcut) 'Rscript.exe'
    if (Test-Path -LiteralPath $candidate -PathType Leaf) { return $candidate }
  }
  $cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
  if (-not $cmd) { $cmd = Get-Command Rscript -ErrorAction SilentlyContinue }
  if ($cmd) { return $cmd.Source }
  throw 'Rscript was not found. Pass -RscriptPath C:\path\to\Rscript.exe.'
}

$Rscript = Resolve-Rscript -ExplicitRscript $RscriptPath -ExplicitRShortcut $RShortcutPath
Write-Host "[singlecell] Skill root: $SkillRoot"
Write-Host "[singlecell] Config:     $Config"
Write-Host "[singlecell] Mode:       $Mode"
Write-Host "[singlecell] Rscript:    $Rscript"

if (-not $SkipDependencyCheck) {
  $depArgs = @((Join-Path $SkillRoot 'install_dependencies.R'))
  if ($InstallPackages) { $depArgs += '--install' }
  if ($InstallGitHubPackages) { $depArgs += '--install-github' }
  if ($RepairCopykat) { $depArgs += '--repair-copykat' }
  & $Rscript @depArgs
  if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

$runArgs = @((Join-Path $SkillRoot 'run_all.R'), $Config, "--mode=$Mode")
$runOutput = & $Rscript @runArgs 2>&1 | Tee-Object -Variable captured
$exitCode = $LASTEXITCODE
$runOutput | ForEach-Object { Write-Host $_ }
if ($exitCode -ne 0) {
  Write-Error "Workflow failed with R exit code $exitCode."
  exit $exitCode
}

$resultLine = $captured | Where-Object { "$_" -match '^RUN_COMPLETE=' } | Select-Object -Last 1
if (-not $resultLine) {
  Write-Error 'R completed without emitting RUN_COMPLETE=<result_root>.'
  exit 20
}
$resultRoot = ("$resultLine" -replace '^RUN_COMPLETE=', '').Trim()
$statusPath = Join-Path $resultRoot '14_qa\FINAL_STATUS.txt'
if (-not (Test-Path -LiteralPath $statusPath -PathType Leaf)) {
  Write-Error "Final QA status file is missing: $statusPath"
  exit 21
}
$status = (Get-Content -LiteralPath $statusPath -Raw).Trim()
if ($status -notin @('PASS','PASS_WITH_WARNINGS')) {
  Write-Error "Final QA did not pass: $status"
  exit 22
}
Write-Host "[singlecell] FINAL_STATUS=$status"
Write-Host "[singlecell] RESULT_ROOT=$resultRoot"
exit 0
