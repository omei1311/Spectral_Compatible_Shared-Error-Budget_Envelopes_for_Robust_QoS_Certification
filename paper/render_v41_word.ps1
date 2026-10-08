$ErrorActionPreference = 'Stop'
$reviewInput = Join-Path $PSScriptRoot 'Spectral_Compatible_Shared_Error_Budget_QoS_Certification_V4.1_EN.docx'
$reviewDir = Join-Path $PSScriptRoot 'results_review/v41_render'
New-Item -ItemType Directory -Path $reviewDir -Force | Out-Null
$reviewPdf = Join-Path $reviewDir 'V4.1_EN.pdf'
$reviewWord = New-Object -ComObject Word.Application
$reviewWord.Visible = $false
$reviewWord.DisplayAlerts = 0
$reviewDoc = $null
try {
    $reviewDoc = $reviewWord.Documents.Open($reviewInput, $false, $true)
    $reviewDoc.Repaginate()
    $reviewDoc.ExportAsFixedFormat($reviewPdf, 17)
    'Pages: ' + $reviewDoc.ComputeStatistics(2)
    'Equations: ' + $reviewDoc.OMaths.Count
    'Tables: ' + $reviewDoc.Tables.Count
} finally {
    if ($null -ne $reviewDoc) { $reviewDoc.Close(0); [Runtime.InteropServices.Marshal]::ReleaseComObject($reviewDoc) | Out-Null }
    if ($reviewWord.Documents.Count -eq 0) { $reviewWord.Quit() }
    [Runtime.InteropServices.Marshal]::ReleaseComObject($reviewWord) | Out-Null
}
& 'C:/Users/57756/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe' -r 120 -png $reviewPdf (Join-Path $reviewDir 'page')
