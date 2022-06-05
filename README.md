# HPC-Final-Project

- You can see the problem description in  doc/Problem.tex
- run code:
```shell
# explicit scheme
make IMPLICIT=no final.out -restart 0
bsub< exec.lsf
# implicit scheme
make IMPLICIT=yes final.out -restart 0
bsub< exec.lsf
```

