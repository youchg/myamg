function A = Generate_dmatcsr_FDM2D_bnd(N, mode)
%N = 4;
%A = Generate_dmatcsr_FDM2D_bnd(5, 4);
%Write_dmatcsr(A, '../dat/
%P = Read_dmatcsr('../output/sum0_64x64_P1.dmatcsr');
%fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
%fd2 = find(~(sum(P(fd1,:),2)==0))
stencil = [-1/2 -2 -1/2; -1 29/4 -1; -1/8 -2 -1/8];
A = Generate_dmatcsr_FDM2D_9P(N, stencil);
switch mode
    case 1
        %=========== case 1: bnd ===============
        for i = 1:N
            A(N*i, :) = zeros(1, N*N);
            A(:, N*i) = zeros(N*N, 1);
            A(N*i, N*i) = 1;
            
            A(N*(i-1)+1, :) = zeros(1, N*N);
            A(:, N*(i-1)+1) = zeros(1, N*N);
            A(N*(i-1)+1, N*(i-1)+1) = 1;
            
            A(i, :) = zeros(1, N*N);
            A(:, i) = zeros(1, N*N);
            A(i, i) = 1;
            
            A(N*(N-1)+i, :) = zeros(1, N*N);
            A(:, N*(N-1)+i) = zeros(1, N*N);
            A(N*(N-1)+i, N*(N-1)+i) = 1;
        end
    case 2
        %=========== case 2: row sum = 0 ===============
        for i = 1:N*N
            s = sum(A(i,:));
            if s==0
                continue;
            else
                aa = A(i, i);
                ab = s-aa;
                A(i, i) = -ab;
            end
        end
    case 3
        %=========== case 3: row sum = 0 with zero row ===============
        for i = 1:N*N
            s = sum(A(i,:));
            if s==0
                continue;
            else
                A(i,:) = zeros(1, N*N);
            end
        end
    case 4
        %=========== case 4: row sum = 0 with zero row and col ===============
        index = [];
        for i = 1:N*N
            s = sum(A(i,:));
            if s==0
                continue;
            else
                index = [index, i];
            end
        end
        A(index, :) = zeros(length(index), N*N);
        A(:, index) = zeros(N*N, length(index));                     
end
