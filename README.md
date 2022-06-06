# HPC-Final-Project

- You can see the problem description in  doc/Problem.tex
- run code:
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


