; GeneScreen 3.0 Inno Setup 脚本
; 用法: 用 Inno Setup Compiler 编译此文件生成 setup.exe

#define MyAppName "GeneScreen"
#define MyAppVersion "3.0"
#define MyAppPublisher "xbzhang"
#define MyAppExeName "GeneScreen.exe"

[Setup]
AppId={{B8F3A2D1-5E7C-4A9B-8D6F-1C2E3A4B5C6D}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
AllowNoIcons=yes
; 输出目录和文件名
OutputDir=dist
OutputBaseFilename=GeneScreen_Setup_{#MyAppVersion}
; 图标（可选）
; SetupIconFile=ui\resources\icons\app.ico
Compression=lzma2
SolidCompression=yes
WizardStyle=modern

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
; 打包 Nuitka 生成的整个目录
Source: "dist\main.dist\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\卸载 {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent
