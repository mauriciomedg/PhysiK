@echo off
setlocal

cd /d "%~dp0"

set UNITY_ROOT=..\PhysiKUnity
set UNITY_PLUGIN_DIR=%UNITY_ROOT%\Assets\Plugins
set CONFIG=Release

if exist build rmdir /s /q build

cmake -S . -B build -DPHYSIK_BUILD_TESTS=ON
cmake --build build --config %CONFIG%

if exist "%UNITY_ROOT%" (
    if not exist "%UNITY_PLUGIN_DIR%" mkdir "%UNITY_PLUGIN_DIR%"

    copy /Y "build\%CONFIG%\Physik.dll" "%UNITY_PLUGIN_DIR%\Physik.dll"

    if exist "build\%CONFIG%\Physik.pdb" (
        copy /Y "build\%CONFIG%\Physik.pdb" "%UNITY_PLUGIN_DIR%\Physik.pdb"
    )
) else (
    echo Unity folder not found: %UNITY_ROOT%
    echo Skipping DLL copy.
)

ctest --test-dir build -C %CONFIG% --output-on-failure

pause