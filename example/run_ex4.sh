#/bin/bash

#--------------------------------------------------------------
for((nev=16;nev<=30;nev++))
do
    ./ex4_10_1_10_1.exe -ini ini/ex4_fem2d_poisson_square_jump_10_1_10_1.dat -nev $nev >result-2016-0406/re_ex4_c8515_nev${nev}.txt 2>&1 &
done
