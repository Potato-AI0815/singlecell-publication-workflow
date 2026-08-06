@echo off
setlocal
set "SKILL_ROOT=%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%SKILL_ROOT%scripts\run_one_click.ps1" %*
set "EXIT_CODE=%ERRORLEVEL%"
endlocal & exit /b %EXIT_CODE%
