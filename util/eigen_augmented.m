function [eval, evec] = eigen_augmented(A, M, P, R, nev, V, tol, ncorrection, eig_val)
% Eigen solver using augmented subspace algorithm.
% Input:
%         A           -- stiff         matrices
%         M           -- mass          matrices
%         P           -- interpolation matrices
%         R           -- restriction   matrices
%         tol         -- tolerance
%         nev         -- number of smallest eigenpairs
%         ncorrection -- max number of iterations
% Output:
%         eval        -- eigenvalues
%         evec        -- eigenvectors
coarsest_level = size(A, 1);
nlevel         = size(A, 1);

Ah = A{1};
Mh = M{1};
AH = A{coarsest_level};
MH = M{coarsest_level};

for iter = 1 : ncorrection
    fprintf(1, '======================== iter %04d ========================\n', iter);
    %===================================
    RAhV = Ah * V;
    for j = 1 : nlevel-1
        RAhV = R{j} * RAhV;
    end
    AHL = [AH RAhV; RAhV' V'*Ah*V];
    
    RMhV = Mh * V;
    for j = 1 : nlevel-1
        RMhV = R{j} * RMhV;
    end
    MHL = [MH RMhV; RMhV' V'*Mh*V];
    
    [VV, D] = eigs(AHL, MHL, nev, 'sm');
    [D, index] = sort(diag(D), 'ascend');
    VV = VV(:, index);
    
    for j = 1 : nev
        rayleigh_quotient(j) = VV(:, j)' * AHL * VV(:, j) / (VV(:, j)' * MHL * VV(:, j));
        fprintf(1, '%04d  %20.15f  %20.15f  %20.15f\n', j, D(j), rayleigh_quotient(j), D(j)-rayleigh_quotient(j));
    end
    fprintf(1, '\n');
    
    %===================================
    PVV = VV(1:end-nev, :);
    for j = nlevel-1 : -1 : 1
        PVV = P{j} * PVV;
    end
    if iter == 1
        V = PVV;
    else
        V = PVV + V * diag(D);
    end
    
    for j = 1 : nev
        rayleigh_quotient(j) = V(:, j)' * Ah * V(:, j) / (V(:, j)' * Mh * V(:, j));
        fprintf(1, '%04d  %20.15f  %20.15f  %20.15f\n', j, D(j), rayleigh_quotient(j), D(j)-rayleigh_quotient(j));
    end
    fprintf(1, '\n');
    
    %===================================
    V = Ah \ (Mh * V * diag(D));
    
    %===================================
    sum_resi = 0;
    rayleigh_quotient = zeros(nev, 1);
    resi              = zeros(nev, 1);
    for j = 1 : nev
        rayleigh_quotient(j) = V(:, j)' * Ah * V(:, j) / (V(:, j)' * Mh * V(:, j));
        %resi(j)              = norm(Ah * V(:,j) - rayleigh_quotient(j) * Mh * V(:, j));
        resi(j)              = norm(Ah * V(:,j) - D(j) * Mh * V(:, j));
        fprintf(1, '%04d  %20.15f  %20.15f  %20.15f  %20.15f\n', j, D(j), rayleigh_quotient(j), resi(j), D(j)-eig_val(j));
        sum_resi = sum_resi + resi(j);
    end
    if resi < tol
        break;
    end
end

eval = rayleigh_quotient;
evec = V;

end