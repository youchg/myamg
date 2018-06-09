clear; 
format long;
%clc;
for i = 1:15
    fprintf(1, '\n');
end

level_finest   = 7;
level_coarsest = 4;
level_number   = level_finest - level_coarsest + 1;

nev = 6;
correction = 5;

file_prefix = '../../dat/tmp/';

disp('======================== read matrices =======================');
disp('Read matrix A...');
A    = cell(level_number, 1);
A{1} = Read_dmatcsr([file_prefix, 'gmg_A_refine7.dat']);
A{2} = Read_dmatcsr([file_prefix, 'gmg_A_refine6.dat']);
A{3} = Read_dmatcsr([file_prefix, 'gmg_A_refine5.dat']);
A{4} = Read_dmatcsr([file_prefix, 'gmg_A_refine4.dat']);
%A{5} = Read_dmatcsr([file_prefix, 'gmg_A_refine3.dat']);
disp('Read matrix A... done.');

disp('Read matrix M...');
M    = cell(level_number, 1);
M{1} = Read_dmatcsr([file_prefix, 'gmg_M_refine7.dat']);
M{2} = Read_dmatcsr([file_prefix, 'gmg_M_refine6.dat']);
M{3} = Read_dmatcsr([file_prefix, 'gmg_M_refine5.dat']);
M{4} = Read_dmatcsr([file_prefix, 'gmg_M_refine4.dat']);
%M{5} = Read_dmatcsr('gmg_M_refine3.dat']);
disp('Read matrix M... done.');

disp('Read matrix P...');
P    = cell(level_number-1, 1);
P{1} = Read_dmatcsr([file_prefix, 'gmg_P_refine7.dat']);
P{2} = Read_dmatcsr([file_prefix, 'gmg_P_refine6.dat']);
P{3} = Read_dmatcsr([file_prefix, 'gmg_P_refine5.dat']);
%P{4} = Read_dmatcsr([file_prefix, 'gmg_P_refine4.dat']);
disp('Read matrix P... done.');

disp('Read matrix R...');
R    = cell(level_number-1, 1);
R{1} = Read_dmatcsr([file_prefix, 'gmg_R_refine7.dat']);
R{2} = Read_dmatcsr([file_prefix, 'gmg_R_refine6.dat']);
R{3} = Read_dmatcsr([file_prefix, 'gmg_R_refine5.dat']);
%R{4} = Read_dmatcsr('[file_prefix, gmg_R_refine4.dat']);
disp('Read matrix R...');
disp('==============================================================');

disp('  ');
dof_finest = size(A{1}, 1);
V = rand(dof_finest, nev);
%eig_val = sort(eigs(A{1}, M{1}, nev, 'sm'), 'ascend');
eig_val = [19.742181652313466
           49.360802349344851
           49.367944185095844
           79.004391701676624
           98.754512911504392
           98.754533209295090];
[eval, evec] = eigen_augmented(A, M, P, R, nev, V, 1.0e-12, 4, eig_val);

% disp('==================== coarsest eigen value ====================');
% AH = A{level_coarsest};
% MH = M{level_coarsest};
% [V, eval] = eigs(AH, MH, nev, 'sm');
% [eval, index] = sort(diag(eval), 'ascend');
% V = V(:, index);
% eval
% disp('==============================================================');
% disp('');

% norm_val = zeros(nev, 2);
% for j = 1:nev
%     norm_val(j,:) = [norm(V(:, j)), norm(P{3} * V(:, j))];
% end
% norm_val

% Ah = A{3};
% Mh = M{3};
% V  = P{3} * V;
% V  = Ah \ (Mh * V);
% norm_val = zeros(nev, 1);
% for j = 1:nev
%     norm_val(j,:) = [norm(V(:, j))];
% end
% norm_val
%
% RAhV = R{3} * Ah * V;
% RMhV = R{3} * Mh * V;
% VTAhV = V'*Ah*V;
% VTMhV = V'*Mh*V;
% RAhV;
% RMhV;
% AHL = [AH, RAhV; RAhV', V'*Ah*V];
% MHL = [MH, RMhV; RMhV', V'*Mh*V];
% [VV, eval] = eigs(AHL, MHL, nev, 'sm');
% [eval, index] = sort(diag(eval), 'ascend');
% V = V(:, index);
% eval

%AHL(:, end-nev+1:end)
%AHL2(:, end-nev+1:end)
