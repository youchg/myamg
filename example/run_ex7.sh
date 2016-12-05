#/bin/bash
echo 'c2458 13'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay2.dat -nb 0 -ne 12 -nev 13 | tee re_ex7_c2458_nev13_nev13.txt
echo 'c2458 14'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay2.dat -nb 0 -ne 12 -nev 14 | tee re_ex7_c2458_nev13_nev14.txt
echo 'c2458 15'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay2.dat -nb 0 -ne 12 -nev 15 | tee re_ex7_c2458_nev13_nev15.txt
echo 'c2458 16'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay2.dat -nb 0 -ne 12 -nev 16 | tee re_ex7_c2458_nev13_nev16.txt
echo 'c2458 17'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay2.dat -nb 0 -ne 12 -nev 17 | tee re_ex7_c2458_nev13_nev17.txt

echo 'c865 13'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay.dat -nb 0 -ne 12 -nev 13 | tee re_ex7_c865_nev13_nev13.txt
echo 'c865 14'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay.dat -nb 0 -ne 12 -nev 14 | tee re_ex7_c865_nev13_nev14.txt
echo 'c865 15'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay.dat -nb 0 -ne 12 -nev 15 | tee re_ex7_c865_nev13_nev15.txt
echo 'c865 16'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay.dat -nb 0 -ne 12 -nev 16 | tee re_ex7_c865_nev13_nev16.txt
echo 'c865 17'
./ex7_lshape_delaunay.exe -ini ini/ex7_fem2d_poisson_lshape_delaunay.dat -nb 0 -ne 12 -nev 17 | tee re_ex7_c865_nev13_nev17.txt


echo 'done'
