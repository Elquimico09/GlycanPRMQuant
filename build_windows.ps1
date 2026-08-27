param(
    [string]$EnvironmentName = "defaultenv"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$projectRoot = (Resolve-Path -LiteralPath $PSScriptRoot).Path
$packagingRoot = [System.IO.Path]::GetFullPath((Join-Path $projectRoot ".packaging"))
$projectPrefix = $projectRoot.TrimEnd([System.IO.Path]::DirectorySeparatorChar) + [System.IO.Path]::DirectorySeparatorChar

if (-not $packagingRoot.StartsWith($projectPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Packaging output resolved outside the project: $packagingRoot"
}

$condaCommand = Get-Command conda -ErrorAction SilentlyContinue
$condaExe = if ($condaCommand -and $condaCommand.CommandType -eq "Application") {
    $condaCommand.Source
} else {
    @(
        (Join-Path $env:USERPROFILE "anaconda3\Scripts\conda.exe")
        (Join-Path $env:USERPROFILE "miniconda3\Scripts\conda.exe")
    ) | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}

if (-not $condaExe) {
    throw "Could not locate conda.exe. Initialize Conda or update build_windows.ps1 with its path."
}

Push-Location $projectRoot
try {
    & $condaExe run -n $EnvironmentName python -c "import PyInstaller"
    if ($LASTEXITCODE -ne 0) {
        throw "PyInstaller is not available in Conda environment '$EnvironmentName'."
    }

    if (Test-Path -LiteralPath $packagingRoot) {
        Remove-Item -LiteralPath $packagingRoot -Recurse -Force
    }

    $buildRoot = New-Item -ItemType Directory -Path (Join-Path $packagingRoot "build") -Force
    $distRoot = New-Item -ItemType Directory -Path (Join-Path $packagingRoot "dist") -Force
    $releaseRoot = New-Item -ItemType Directory -Path (Join-Path $packagingRoot "release") -Force

    & $condaExe run -n $EnvironmentName python -m PyInstaller `
        --noconfirm `
        --clean `
        --workpath $buildRoot.FullName `
        --distpath $distRoot.FullName `
        (Join-Path $projectRoot "GlycanPRMQuant.spec")
    if ($LASTEXITCODE -ne 0) {
        throw "PyInstaller build failed with exit code $LASTEXITCODE."
    }

    $appRoot = Join-Path $distRoot.FullName "GlycanPRMQuant"
    $exePath = Join-Path $appRoot "GlycanPRMQuant.exe"
    if (-not (Test-Path -LiteralPath $exePath)) {
        throw "Expected executable was not produced: $exePath"
    }

    $version = ((& $condaExe run -n $EnvironmentName python -c "import tomllib; print(tomllib.load(open('pyproject.toml', 'rb'))['project']['version'], end='')") -join "").Trim()
    if ($LASTEXITCODE -ne 0 -or -not $version) {
        throw "Could not read the project version from pyproject.toml."
    }

    $archivePath = Join-Path $releaseRoot.FullName "GlycanPRMQuant-Windows-x64-v$version.zip"
    Compress-Archive -Path (Join-Path $appRoot "*") -DestinationPath $archivePath -CompressionLevel Optimal

    $hash = Get-FileHash -LiteralPath $archivePath -Algorithm SHA256
    $hashLine = "$($hash.Hash)  $([System.IO.Path]::GetFileName($archivePath))"
    Set-Content -LiteralPath "$archivePath.sha256" -Value $hashLine -Encoding ascii

    Write-Host "Build complete:"
    Write-Host "  Application: $appRoot"
    Write-Host "  Release ZIP: $archivePath"
    Write-Host "  SHA-256:     $($hash.Hash)"
    Write-Host "Scan the ZIP and extracted application with Microsoft Defender before publishing."
} finally {
    Pop-Location
}
