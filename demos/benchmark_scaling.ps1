<#
.SYNOPSIS
  Parallel scaling benchmark for e-XSpice: dense AC sweep on an RC ladder.

.DESCRIPTION
  One operating point and one linearization, then 5001 complex LU solves.
  The frequency loop is embarrassingly parallel, so this is the workload
  that shows the OpenMP (threads) and MPI (.OPTIONS ACFREQPAR) paths.

  Records the simulation time reported by the simulator and the SHA-256 of
  every CSV, so the byte-identity of the distributed runs can be checked
  from the artifacts themselves. Timings are informative, not contractual.

.PARAMETER Exe
  Path to the e-XSpice executable (default ../build/e-XSpice.exe).

.PARAMETER Csv2Svg
  Path to csv2svg; when present, also rebuilds demos/plots/speedup.svg.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File demos/benchmark_scaling.ps1
#>

param(
    [string]$Exe = "$PSScriptRoot/../build/e-XSpice.exe",
    [string]$Csv2Svg = "$PSScriptRoot/../build/csv2svg.exe",
    [string]$OutDir = "$PSScriptRoot/out"
)

$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force $OutDir | Out-Null

$circuit = "$PSScriptRoot/ac_ladder_200_freqpar.cir"

$configs = @(
    @{ Name = '1 thread';    Omp = 1; Ranks = 1 },
    @{ Name = '2 threads';   Omp = 2; Ranks = 1 },
    @{ Name = '4 threads';   Omp = 4; Ranks = 1 },
    @{ Name = '8 threads';   Omp = 8; Ranks = 1 },
    @{ Name = 'MPI 2 ranks'; Omp = 1; Ranks = 2 },
    @{ Name = 'MPI 4 ranks'; Omp = 1; Ranks = 4 }
)

$rows = @()
foreach ($cfg in $configs) {
    $tag = $cfg.Name -replace '\s+', ''
    $csv = Join-Path $OutDir "scaling_$tag.csv"
    $env:OMP_NUM_THREADS = [string]$cfg.Omp

    $sw = [Diagnostics.Stopwatch]::StartNew()
    if ($cfg.Ranks -gt 1) {
        $log = & mpiexec -n $cfg.Ranks $Exe file $circuit $csv 2>&1 | Out-String
    } else {
        $log = & $Exe file $circuit $csv 2>&1 | Out-String
    }
    $sw.Stop()
    $rc = $LASTEXITCODE

    $m = [regex]::Match($log, '\|\s*time:\s*([\d.eE+-]+)')
    $time = if ($m.Success) { [double]$m.Groups[1].Value } else { [double]::NaN }
    $hash = if (Test-Path $csv) { (Get-FileHash $csv -Algorithm SHA256).Hash } else { '-' }

    $rows += [pscustomobject]@{
        Config  = $cfg.Name
        Threads = $cfg.Omp
        Ranks   = $cfg.Ranks
        Time_s  = $time
        Wall_s  = [math]::Round($sw.Elapsed.TotalSeconds, 3)
        Hash    = $hash
        Exit    = $rc
    }
    "{0,-14} exit={1} time={2,8:N3}s wall={3,8:N3}s hash={4}" -f `
        $cfg.Name, $rc, $time, $sw.Elapsed.TotalSeconds, $hash.Substring(0, 12)

    # Let the CPU settle between configurations: laptop boost/thermal state
    # otherwise dominates the differences this script is trying to show.
    Start-Sleep -Seconds 4
}

# Write CSVs by hand with INVARIANT number formatting: Export-Csv follows the
# machine locale (a pt-BR Windows would emit decimal commas) and csv2svg reads
# dot-decimal fields.
$inv = [System.Globalization.CultureInfo]::InvariantCulture
$resLines = @('Config,Threads,Ranks,Time_s,Wall_s,Hash,Exit')
foreach ($r in $rows) {
    $resLines += '{0},{1},{2},{3},{4},{5},{6}' -f `
        $r.Config, $r.Threads, $r.Ranks, `
        $r.Time_s.ToString('F3', $inv), $r.Wall_s.ToString('F3', $inv), `
        $r.Hash, $r.Exit
}
Set-Content -Path (Join-Path $OutDir 'scaling_results.csv') -Value $resLines

$base = ($rows | Where-Object { $_.Config -eq '1 thread' }).Time_s
$speedRows = foreach ($r in $rows) {
    [pscustomobject]@{
        config  = $r.Config
        time_s  = $r.Time_s
        speedup = [math]::Round($base / $r.Time_s, 2)
    }
}
$speedLines = @('config,time_s,speedup')
foreach ($r in $speedRows) {
    $speedLines += '{0},{1},{2}' -f `
        $r.config, $r.time_s.ToString('F3', $inv), $r.speedup.ToString('F2', $inv)
}
Set-Content -Path (Join-Path $OutDir 'scaling_speedup.csv') -Value $speedLines

Write-Host ''
Write-Host 'Speedups (relative to 1 thread):'
foreach ($r in $speedRows) { "  {0,-14} {1,6:N2}x   ({2:N3}s)" -f $r.config, $r.speedup, $r.time_s }

$serialHash = ($rows | Where-Object { $_.Config -eq '1 thread' }).Hash
$mpiHashes  = $rows | Where-Object { $_.Ranks -gt 1 } | ForEach-Object { $_.Hash }
if (($mpiHashes | Where-Object { $_ -ne $serialHash }).Count -eq 0) {
    Write-Host ''
    Write-Host 'CSV byte-identical between serial and every MPI rank count.'
} else {
    Write-Host ''
    Write-Host 'WARNING: some MPI run produced a different CSV (see hashes above).'
}

if ($Csv2Svg -and (Test-Path $Csv2Svg)) {
    $plot = "$PSScriptRoot/plots/speedup.svg"
    New-Item -ItemType Directory -Force (Split-Path $plot) | Out-Null
    & $Csv2Svg (Join-Path $OutDir 'scaling_speedup.csv') $plot `
        --type bar --series 2 `
        --title "AC sweep scaling - 5001 frequencies, 202 unknowns" `
        --ylabel "Speedup (x)"
    Write-Host "figure: $plot"
}
