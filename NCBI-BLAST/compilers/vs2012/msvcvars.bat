@echo off
REM
REM $Id: msvcvars.bat 407440 2013-07-22 13:00:11Z gouriano $
REM

@if not "%VSINSTALLDIR%"=="" goto devenv
@call "%VS110COMNTOOLS%vsvars32.bat"

:devenv

if exist "%VS110COMNTOOLS%..\IDE\VCExpress.*" set DEVENV="%VS110COMNTOOLS%..\IDE\VCExpress"
if exist "%VS110COMNTOOLS%..\IDE\devenv.*" set DEVENV="%VS110COMNTOOLS%..\IDE\devenv"

:end
