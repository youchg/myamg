#/bin/bash

#--------------------------------------------------------------
#for((nev=1;nev<=9;nev++))
#do
#    ./ex2_lshape.exe -ini ini/ex2_fem2d_poisson_lshape.dat -nev $nev >result-2016-0406/re_ex2_nev0${nev}_c500.txt 2>&1 &
#done

for((nev=16;nev<=30;nev++))
do
    ./ex2_lshape.exe -ini ini/ex2_fem2d_poisson_lshape.dat -nev $nev >result-2016-0406/re_ex2_nev${nev}_c413.txt 2>&1 &
done

