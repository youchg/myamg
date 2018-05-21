clear; clc;

level_finest   = 7;
level_coarsest = 4;
level_number   = level_finest - level_coarsest + 1;

nev = 5;
correction = 5;

A    = cell(level_number, 1);
A{1} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_A_refine7.dat');
A{2} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_A_refine6.dat');
A{3} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_A_refine5.dat');
A{4} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_A_refine4.dat');
%A{5} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_A_refine3.dat');

M    = cell(level_number, 1);
M{1} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_M_refine7.dat');
M{2} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_M_refine6.dat');
M{3} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_M_refine5.dat');
M{4} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_M_refine4.dat');
%M{5} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_M_refine3.dat');

P    = cell(level_number-1, 1);
P{1} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_P_refine7.dat');
P{2} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_P_refine6.dat');
P{3} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_P_refine5.dat');
%P{4} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_P_refine4.dat');

R    = cell(level_number-1, 1);
R{1} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_R_refine7.dat');
R{2} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_R_refine6.dat');
R{3} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_R_refine5.dat');
%R{4} = Read_dmatcsr('/Users/ycg/Mytest/tmp/gmg_R_refine4.dat');

disp('==================== coarsest eigen value ====================');
AH = A{4};
MH = M{4};
[V, eval] = eigs(AH, MH, nev, 'sm');
[eval, index] = sort(diag(eval), 'ascend');
V = V(:, index);
eval

% norm_val = zeros(nev, 2);
% for j = 1:nev
%     norm_val(j,:) = [norm(V(:, j)), norm(P{3} * V(:, j))];
% end
% norm_val

Ah = A{3};
Mh = M{3};
V  = P{3} * V;
V  = Ah \ (Mh * V);
norm_val = zeros(nev, 1);
for j = 1:nev
    norm_val(j,:) = [norm(V(:, j))];
end
norm_val

RAhV = R{3} * Ah * V;
RMhV = R{3} * Mh * V;
VTAhV = V'*Ah*V
VTMhV = V'*Mh*V
RAhV
RMhV
AHL = [AH, RAhV; RAhV', V'*Ah*V];
MHL = [MH, RMhV; RMhV', V'*Mh*V];
[VV, eval] = eigs(AHL, MHL, nev, 'sm');
[eval, index] = sort(diag(eval), 'ascend');
V = V(:, index);
eval

AHL2 = Read_dmatcsr('/Users/ycg/Mytest/Alarge_seq.dat');
MHL2 = Read_dmatcsr('/Users/ycg/Mytest/Mlarge_seq.dat');
norm(AHL-AHL2)
norm(MHL-MHL2)

%AHL(:, end-nev+1:end)
%AHL2(:, end-nev+1:end)


% Ah = A{1};
% Mh = M{1};
% AH = A{level_number};
% MH = M{level_number};

%dof_finest = size(A{1}, 1);
%V = rand(dof_finest, nev);

% for i = 1 : 1
%     AH   = A{level_number};
%     RAhV = Ah * V;
%     for j = 1 : level_number-1
%         RAhV = R{j} * RAhV;
%     end
%     AHL = [AH RAhV; RAhV' V'*Ah*V];
%     
%     MH   = M{level_number};
%     RMhV = Mh * V;
%     for j = 1 : level_number-1
%         RMhV = R{j} * RMhV;
%     end
%     MHL = [MH RMhV; RMhV' V'*Mh*V];
%     
%     [evec, eval] = eigs(AHL, MHL, nev, 'sm');
%     eval = diag(eval)
% end