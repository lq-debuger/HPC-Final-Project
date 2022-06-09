# HPC-Final-Project
- SID:12132402
- NAME：李强
## Prerequisites
- petsc
- gnuplot
- latex

## Compile
```shell
#$ explicit scheme
make IMPLICIT=no final.out -restart 0
## implicit scheme
make IMPLICIT=yes final.out -restart 0
```
- If you want to start with the existing H5 file, add -restart 1 when compiling, otherwise add -restart 0

## Run
```shell
bsub< exec.lsf
```

## Document 
- The problem description in  doc/Problem.tex
- The report is in report/report.tex
- Imgs and relative data are in plot
- The calculting result is in log


