#/bin/bash

for i in 2 5 7 9
do
    ./ex3_1000_1_001_1.exe -ini ini/ex3_fem2d_poisson_square_jump_1000_1_001_1.dat -nev $i >result-2016-0406/re_ex3_nev0$i.txt 2>&1 &
done

for i in 12 14 16 18 21 23 25 27 29
do
    ./ex3_1000_1_001_1.exe -ini ini/ex3_fem2d_poisson_square_jump_1000_1_001_1.dat -nev $i >result-2016-0406/re_ex3_nev$i.txt 2>&1 &
done

