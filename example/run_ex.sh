#/bin/bash

#--------------------------------------------------------------
#for((nev=1;nev<=9;nev++))
#do
#    ./ex2_lshape.exe -ini ini/ex2_fem2d_poisson_lshape.dat -nev $nev >result-2016-0406/re_ex2_nev0${nev}_c500.txt 2>&1 &
#done
#
#for((nev=10;nev<=30;nev++))
#do
#    ./ex2_lshape.exe -ini ini/ex2_fem2d_poisson_lshape.dat -nev $nev >result-2016-0406/re_ex2_nev${nev}_c500.txt 2>&1 &
#done

#--------------------------------------------------------------
for((nev=1;nev<=9;nev++))
do
    ./ex3_1000_1_001_1.exe -ini ini/ex3_fem2d_poisson_square_jump_1000_1_001_1.dat -nev $nev >result-2016-0406/re_ex3_nev0${nev}_c500.txt 2>&1 &
done

for((nev=10;nev<=30;nev++))
do
    ./ex3_1000_1_001_1.exe -ini ini/ex3_fem2d_poisson_square_jump_1000_1_001_1.dat -nev $nev >result-2016-0406/re_ex3_nev${nev}_c500.txt 2>&1 &
done

#--------------------------------------------------------------
for((nev=1;nev<=9;nev++))
do
    ./ex3_1000_1_001_1.exe -ini ini/ex3_fem2d_poisson_square_jump_1000_1_001_1.dat -nev $nev >result-2016-0406/re_ex4_nev0${nev}_c500.txt 2>&1 &
done

for((nev=10;nev<=30;nev++))
do
    ./ex3_1000_1_001_1.exe -ini ini/ex3_fem2d_poisson_square_jump_1000_1_001_1.dat -nev $nev >result-2016-0406/re_ex4_nev${nev}_c500.txt 2>&1 &
done
