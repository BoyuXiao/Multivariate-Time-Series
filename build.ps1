# 一键编译 TLB（Visual Studio 2022）
$ErrorActionPreference = "Stop"
$ProjectRoot = $PSScriptRoot
$BuildDir = Join-Path $ProjectRoot "build"

if (-not (Test-Path $BuildDir)) {
    New-Item -ItemType Directory -Path $BuildDir | Out-Null
}

Set-Location $BuildDir

if (-not (Test-Path "CMakeCache.txt")) {
    Write-Host ">>> CMake configure ..."
    cmake .. -G "Visual Studio 17 2022" -A x64
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

Write-Host ">>> Building Release ..."
cmake --build . --config Release
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

$Exe = Join-Path $BuildDir "Release\TLB.exe"
if (Test-Path $Exe) {
    Write-Host ">>> OK: $Exe"
} else {
    Write-Host ">>> Build finished but TLB.exe not found under Release/"
    exit 1
}
