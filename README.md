# HPC-Final-Project

## Compile
```shell
#$ explicit scheme
make IMPLICIT=no final.out -restart 0
bsub< exec.lsf
## implicit scheme
make IMPLICIT=yes final.out -restart 0
bsub< exec.lsf
## build options
```
- If you want to start with the existing H5 file, add -restart 1 when compiling, otherwise add -restart 0

## Document 
- The problem description in  doc/Problem.tex
- The report is in report/report.tex
- Imgs and relative data are in plot
- The calculting result is in log


