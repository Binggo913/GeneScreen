param(
    [string]$SourceDir,
    [string]$OutputDir,
    [string]$AppVersion = "1.0",
    [string]$Arch = "x64",
    [string]$Timestamp,
    [string]$IsccPath
)

$ErrorActionPreference = "Stop"

$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$InstallerScript = Join-Path $ProjectDir "installer\GeneScreen.iss"

if (-not $SourceDir) {
    $SourceDir = Get-ChildItem -Path (Join-Path $ProjectDir "dist") -Directory -Filter "*.dist" |
        Select-Object -First 1 -ExpandProperty FullName
}
if (-not $SourceDir -or -not (Test-Path $SourceDir)) {
    throw "Nuitka .dist source directory not found. Run build.py first, or pass -SourceDir."
}

if (-not $OutputDir) {
    $OutputDir = Join-Path $ProjectDir "dist\packages"
}
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

if (-not $Timestamp) {
    $Timestamp = Get-Date -Format "yyyyMMdd-HHmmss"
}

function Find-Iscc {
    param([string]$ExplicitPath)

    if ($ExplicitPath) {
        if (Test-Path $ExplicitPath) {
            $resolvedPath = (Resolve-Path $ExplicitPath).Path
            if ([System.IO.Path]::GetExtension($resolvedPath).Equals(".lnk", [System.StringComparison]::OrdinalIgnoreCase)) {
                $shortcut = (New-Object -ComObject WScript.Shell).CreateShortcut($resolvedPath)
                if ($shortcut.TargetPath -and (Test-Path $shortcut.TargetPath)) {
                    if ([System.IO.Path]::GetFileName($shortcut.TargetPath).Equals("ISCC.exe", [System.StringComparison]::OrdinalIgnoreCase)) {
                        return $shortcut.TargetPath
                    }

                    $targetDir = Split-Path -Parent $shortcut.TargetPath
                    $siblingCompiler = Join-Path $targetDir "ISCC.exe"
                    if (Test-Path $siblingCompiler) {
                        return $siblingCompiler
                    }
                }

                throw "ISCC.exe could not be resolved from shortcut: $ExplicitPath"
            }

            return $resolvedPath
        }
        throw "ISCC.exe not found at explicit path: $ExplicitPath"
    }

    $command = Get-Command "ISCC.exe" -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }

    $candidates = @(
        "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
        "${env:ProgramFiles}\Inno Setup 6\ISCC.exe",
        "${env:LOCALAPPDATA}\Programs\Inno Setup 6\ISCC.exe",
        "${env:LOCALAPPDATA}\Programs\Inno Setup\ISCC.exe",
        "D:\Inno Setup 6\ISCC.exe"
    )

    foreach ($candidate in $candidates) {
        if ($candidate -and (Test-Path $candidate)) {
            return $candidate
        }
    }

    throw "ISCC.exe not found. Install Inno Setup 6 or pass -IsccPath."
}

$Iscc = Find-Iscc -ExplicitPath $IsccPath

Write-Host "Using Inno Setup compiler: $Iscc"
Write-Host "Source directory: $SourceDir"
Write-Host "Output directory: $OutputDir"

& $Iscc `
    "/DSourceDir=$SourceDir" `
    "/DOutputDir=$OutputDir" `
    "/DAppVersion=$AppVersion" `
    "/DArch=$Arch" `
    "/DTimestamp=$Timestamp" `
    $InstallerScript

if ($LASTEXITCODE -ne 0) {
    throw "Inno Setup failed with exit code $LASTEXITCODE"
}

$setupName = "GeneScreen_${AppVersion}_windows_${Arch}_setup_${Timestamp}.exe"
$setupPath = Join-Path $OutputDir $setupName
if (-not (Test-Path $setupPath)) {
    throw "Expected setup package was not created: $setupPath"
}

Write-Host "Setup package created: $setupPath"
