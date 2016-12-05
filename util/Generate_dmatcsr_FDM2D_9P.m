function A = Generate_dmatcsr_FDM2D_9P(N, stencil)
%N: inner grid points for each axis
%A = [ B CR
%     CL B CR
%        ......
%           CL  B CR
%              CL  B];  N^2 x N^2
%stencil = [0 -1 0; -1 4 -1; 0 -1 0];
%stencil = [-1/2 -2 -1/2; -1 29/4 -1; -1/8 -2 -1/8];
%A = Generate_dmatcsr_FDM2D_9P(4,stencil);

      
B  = diag(ones(1, N))*stencil(2,2) + diag(ones(1, N-1),1)*stencil(2,1) ...
     + diag(ones(1, N-1), -1)*stencil(2,3);
     
CL = diag(ones(1, N))*stencil(3,2) + diag(ones(1, N-1),1)*stencil(3,1) ...
     + diag(ones(1, N-1), -1)*stencil(3,3);

CR = diag(ones(1, N))*stencil(1,2) + diag(ones(1, N-1),1)*stencil(1,1) ...
     + diag(ones(1, N-1), -1)*stencil(1,3);

H0 = [B CR];
H1 = [CL B CR];
H2 = [CL B];

ZZ = zeros(N, (N-2)*N);
A = [H0, ZZ];
for i = 2:N-1
    ZZ1 = zeros(N, (i-2)*N);
    ZZ2 = zeros(N, (N-i-1)*N);
    A = [A; ZZ1 H1 ZZ2];
end
A = [A; ZZ H2];

    
