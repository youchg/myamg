A0 = Read_dmatcsr('../dat/fem_poisson_refine_5_A.dmatcsr');
A1 = Read_dmatcsr('../output/fem_poisson_refine_5_read1.dmatcsr');
A2 = Read_dmatcsr('../output/fem_poisson_refine_5_read2.dmatcsr');
A3 = Read_dmatcsr('../output/fem_poisson_refine_5_copy1.dmatcsr');
T1 = A0';
T2 = Read_dmatcsr('../output/fem_poisson_refine_5_trans.dmatcsr');

MA = Read_dmatcsr('../dat/fdm2d_9P_bnd_256x256_A.dmatcsr');
MB = Read_dmatcsr('../dat/sum0_256x256_A.dmatcsr');
C1 = MA*MB;
C2 = Read_dmatcsr('../output/multi_256.dmatcsr');

norm1 = norm(A0-A1)
norm2 = norm(A1-A2)
norm3 = norm(A2-A3)
norm4 = norm(A3-A0)
norm5 = norm(T1-T2)
norm6 = norm(C1-C2)
