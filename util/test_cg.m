%A   = Read_dmatcsr('../../../FEM_SOFT/dat/matrix_data_discontinuous_1000_1_001_1/gmg_A_refine5.m');
%A   = Read_dmatcsr('../dat/test_cg/tridiag_2m1m1_128x128.dmatcsr');
%N = 256;
%A = diag(2*ones(N,1))+diag(-1*ones(N-1,1),1)+diag(-1*ones(N-1,1),-1);
%Write_dmatcsr(A, '../dat/test_cg/tridiag_2m1m1_256x256.dmatcsr');
%Write_dmatcsr_for_fasp(A, '../dat/test_cg/discontinuous_1000_1_001_1.fasp');
%Write_dmatcsr(A, '../dat/test_cg/discontinuous_1000_1_001_1.dmatcsr');
A   = Read_dmatcsr('../dat/xufei/stiff_matrix-10.dat');
sol = ones(size(A,1), 1);
b   = A*sol;
x0  = zeros(size(A,1), 1);
x   = cg(A, b, x0, 10^-14, 20, sol);
rn = norm(b-A*x)
er = norm(x-sol)
x1 = A\b;
er2 = norm(sol-x1)
%condA = cond(A)
%ratio = (sqrt(condA)-1)/(sqrt(condA)+1)
%fid = fopen('b.dat', 'w');
%fprintf(fid, '%e\n', b');
%fclose(fid);