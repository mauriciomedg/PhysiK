@echo off
setlocal

cd /d "%~dp0"

set UNITY_PLUGIN_DIR=..\PhysiKUnity\Assets\Plugins

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" --force

if exist build rmdir /s /q build

cmake -S . -B build -DPHYSIK_BUILD_TESTS=ON
cmake --build build --config Release

if not exist "%UNITY_PLUGIN_DIR%" mkdir "%UNITY_PLUGIN_DIR%"

copy /Y "build\Release\Physik.dll" "%UNITY_PLUGIN_DIR%\Physik.dll"

if exist "build\Release\Physik.pdb" (
    copy /Y "build\Release\Physik.pdb" "%UNITY_PLUGIN_DIR%\Physik.pdb"
)

ctest --test-dir build -C Release --output-on-failure

pause