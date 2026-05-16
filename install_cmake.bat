@echo off
setlocal

echo Installing CMake using winget...

winget install Kitware.CMake --accept-package-agreements --accept-source-agreements

echo.
echo Done.
echo Close this terminal and reopen it, then test:
echo cmake --version
echo ctest --version

pause