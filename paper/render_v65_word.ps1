$ErrorActionPreference = 'Stop'
$taskProject = Split-Path $PSScriptRoot -Parent
$taskInput = Join-Path $PSScriptRoot 'Sharp_Robust_QoS_Certification_from_Six-Scalar_Spectral_Summaries_V6.5_EN.docx'
$taskOutput = Join-Path $taskProject 'results_v65/render'
New-Item -ItemType Directory -Path $taskOutput -Force | Out-Null
$taskPdf = Join-Path $taskOutput 'V6.5_EN.pdf'
$taskWord = New-Object -ComObject Word.Application
$taskWord.Visible = $false
$taskWord.DisplayAlerts = 0
$taskDoc = $null
try {
    $taskDoc = $taskWord.Documents.Open($taskInput, $false, $true)
    $taskDoc.Repaginate()
    $taskDoc.ExportAsFixedFormat($taskPdf, 17)
    'Pages: ' + $taskDoc.ComputeStatistics(2)
    'Equations: ' + $taskDoc.OMaths.Count
    'Tables: ' + $taskDoc.Tables.Count
} finally {
    if ($null -ne $taskDoc) { $taskDoc.Close(0); [Runtime.InteropServices.Marshal]::ReleaseComObject($taskDoc) | Out-Null }
    if ($taskWord.Documents.Count -eq 0) { $taskWord.Quit() }
    [Runtime.InteropServices.Marshal]::ReleaseComObject($taskWord) | Out-Null
}
& 'C:/Users/57756/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe' -r 110 -png $taskPdf (Join-Path $taskOutput 'page')
if ($LASTEXITCODE -ne 0) { throw 'Page rasterization failed' }
