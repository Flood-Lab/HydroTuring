@echo off
rem Dev shim for Windows shells: run the CLI without installing the package.
rem The bash shim `ht` does the same on POSIX.
setlocal
set "PYTHONPATH=%~dp0src;%PYTHONPATH%"
python -m hydroturing %*
endlocal
