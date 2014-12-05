@echo off
REM
REM $Id: msvcvars.bat 430635 2014-03-27 17:34:42Z gouriano $
REM

@if not "%VSINSTALLDIR%"=="" goto devenv
@call "%VS120COMNTOOLS%vsvars32.bat"

:devenv

if exist "%VS120COMNTOOLS%..\IDE\VCExpress.*" set DEVENV="%VS120COMNTOOLS%..\IDE\VCExpress"
if exist "%VS120COMNTOOLS%..\IDE\devenv.*" set DEVENV="%VS120COMNTOOLS%..\IDE\devenv"

:end
