#/bin/bash

for i in 5 7 9
do
    ./ex1_poisson.exe -ini ini/ex1_fem2d_poisson_square.dat -nev $i >result-2016-0406/re_ex1_nev0$i.txt 2>&1 &
done

for i in 12 14 16 18 21 23 25 27 29
do
    ./ex1_poisson.exe -ini ini/ex1_fem2d_poisson_square.dat -nev $i >result-2016-0406/re_ex1_nev$i.txt 2>&1 &
done

