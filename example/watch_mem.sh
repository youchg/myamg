#!/bin/bash
MyName=$USER
TotalMem=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')

./ex1_poisson.exe -ini ini/ex1_fem2d_poisson_square.dat -nev 13 -nb 0 -ne 12 >../output/ex1_poisson_amg_refine11.txt &

MyCommand=ex1_poisson.exe #程序名
MyStep=1 #秒
MyTime=0 #秒
MyInfo=$(ps -u $MyName -o comm,time,etime,pcpu,pmem | grep $MyCommand)
while [[ -n $MyInfo ]]
do
    MyMemPercent=$(echo $MyInfo | awk '{print $5}')
    MyCPUPercent=$(echo $MyInfo | awk '{print $4}')
    MyMem=$(echo "scale=2; $TotalMem*$MyMemPercent/100.0/1000.0" | bc)
    MyTime=$(echo "$MyTime+1" | bc)
    echo $MyTime $MyMem
    #echo $MyTime $MyMem $MyCPUPercent
    sleep $MyStep
    MyInfo=$(ps -u $MyName -o comm,time,etime,pcpu,pmem | grep $MyCommand)
done

