#define AppName "GeneScreen"

#ifndef AppVersion
#define AppVersion "1.0"
#endif

#ifndef Arch
#define Arch "x64"
#endif

#ifndef Timestamp
#define Timestamp GetDateTimeString('yyyymmdd-hhnnss', '', '')
#endif

#ifndef SourceDir
#define SourceDir "..\dist\main.dist"
#endif

#ifndef OutputDir
#define OutputDir "..\dist\packages"
#endif

#define OutputBaseName AppName + "_" + AppVersion + "_windows_" + Arch + "_setup_" + Timestamp

[Setup]
AppId={{ACB05B1A-B9F0-4FB7-8F4D-7B349D557E1E}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher=GeneScreen
AppPublisherURL=https://github.com/Binggo913/GeneScreen
AppSupportURL=https://github.com/Binggo913/GeneScreen
AppUpdatesURL=https://github.com/Binggo913/GeneScreen
DefaultDirName={localappdata}\Programs\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
OutputDir={#OutputDir}
OutputBaseFilename={#OutputBaseName}
SetupIconFile=..\ui\resources\icons\app.ico
UninstallDisplayIcon={app}\GeneScreen.exe
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\GeneScreen.exe"; WorkingDir: "{app}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\GeneScreen.exe"; WorkingDir: "{app}"; Tasks: desktopicon

[Run]
Filename: "{app}\GeneScreen.exe"; Description: "{cm:LaunchProgram,{#AppName}}"; Flags: nowait postinstall skipifsilent
