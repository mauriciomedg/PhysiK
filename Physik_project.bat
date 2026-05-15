@echo off
setlocal

cd /d E:\PhysiK

echo ============================
echo Init Intel oneAPI environment
echo ============================

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" --force

echo ============================
echo Clean old build
echo ============================

if exist build rmdir /s /q build

echo ============================
echo Configure build
echo ============================

cmake -S . -B build -DPHYSIK_BUILD_TESTS=ON

echo ============================
echo Build Release
echo ============================

cmake --build build --config Release

echo ============================
echo Test Release
echo ============================

ctest --test-dir build -C Release --output-on-failure

echo ============================
echo Done
echo ============================

pause