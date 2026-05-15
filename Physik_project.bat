@echo off
setlocal

cd /d E:\PhysiK

echo ============================
echo Commit and push current branch
echo ============================

git status
git add .
git commit -m "Add MKL tiny tests"
git push

echo ============================
echo Init Intel oneAPI environment
echo ============================

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" --force

echo ============================
echo Clean old builds
echo ============================

if exist build rmdir /s /q build
if exist build-mkl rmdir /s /q build-mkl

echo ============================
echo Configure normal build
echo ============================

cmake -S . -B build -DPHYSIK_BUILD_TESTS=ON

echo ============================
echo Build normal
echo ============================

cmake --build build --config Release

echo ============================
echo Test normal
echo ============================

ctest --test-dir build -C Release --output-on-failure

echo ============================
echo Done
echo ============================

pause