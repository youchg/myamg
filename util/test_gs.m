%A   = Read_dmatcsr('../dat/xufei/stiff_matrix-5.dat');
A   = Read_dmatcsr('../../../FEM_SOFT/dat/matrix_data_discontinuous_1000_1_001_1/gmg_A_refine6.m');
%A   = sparse(A);
sol = ones(size(A,1), 1);
b   = A*sol;
% x0  = zeros(size(A,1), 1);
% x   = x0;
% norm(b-A*x)
% 
% D = diag(diag(A));
% L = -tril(A, -1);
% U = -triu(A,  1);
% 
% max_iter = 20;
% nr = zeros(max_iter, 1);
% for i = 1 : max_iter
%     x = (D-L)\(U*x + b);
%     nr(i) = norm(b-A*x);
% end
% ratio = zeros(max_iter, 1);
% ratio(2:end) = nr(2:end)./nr(1:end-1);
% format long
% [nr, ratio]

% G = inv(D-L) * U;
% conv_ratio = max(abs(eig(G)))
% -log(conv_ratio)

Write_dmatcsr_for_fasp(A, ...
    '../output/discontinuous_1000_1_001_1_refine6.fasp');
fid = fopen('../output/b.dat', 'w');
fprintf(fid, '%e\n', b');
fclose(fid);