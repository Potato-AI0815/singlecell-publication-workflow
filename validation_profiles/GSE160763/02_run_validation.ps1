param([switch]$InstallPackages,[ValidateSet("smoke","sct")][string]$Mode="smoke")
$ErrorActionPreference="Stop"
$ProjectRoot="C:\Users\YHN\Desktop\Qoder\GSE160763_analysis"
$RawTar="C:\Users\YHN\Desktop\Qoder\GSE160763_RAW.tar"
$Rscript="D:\Ruanjian\R-4.5.3\bin\x64\Rscript.exe"
$Profile=Split-Path -Parent $MyInvocation.MyCommand.Path
$SkillRoot=Split-Path -Parent (Split-Path -Parent $Profile)
if(!(Test-Path $Rscript)){throw "Rscript.exe not found: $Rscript"}
& $Rscript (Join-Path $Profile "00_prepare_GSE160763.R") $ProjectRoot $RawTar
$installArgs=@((Join-Path $Profile "01_install_core_dependencies.R"));if($InstallPackages){$installArgs += "--install"}
& $Rscript @installArgs
$config=if($Mode -eq "smoke"){"config_01_core_smoke.yml"}else{"config_02_core_sct_gate_test.yml"}
& $Rscript (Join-Path $SkillRoot "run_all.R") (Join-Path $Profile $config)
$latest=Get-Content (Join-Path $ProjectRoot "results\LATEST_RUN.txt") -First 1
& $Rscript (Join-Path $Profile "03_validate_GSE160763_outputs.R") $latest $Mode
Write-Host "Validation completed. Result root: $latest"
