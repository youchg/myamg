cc%=========================================================
P = Read_dmatcsr('../output/sum0_16x16_P1.dmatcsr');
fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
fd2 = find(~(sum(P(fd1,:),2)==0))

P = Read_dmatcsr('../output/sum0_64x64_P1.dmatcsr');
fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
fd2 = find(~(sum(P(fd1,:),2)==0))

P = Read_dmatcsr('../output/sum0_256x256_P1.dmatcsr');
fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
fd2 = find(~(sum(P(fd1,:),2)==0))

%=========================================================
P = Read_dmatcsr('../output/sum0_zero_row_16x16_P1.dmatcsr');
fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
fd2 = find(~(sum(P(fd1,:),2)==0))

P = Read_dmatcsr('../output/sum0_zero_row_25x25_P1.dmatcsr');
fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
fd2 = find(~(sum(P(fd1,:),2)==0))

P = Read_dmatcsr('../output/sum0_zero_row_64x64_P1.dmatcsr');
fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
fd2 = find(~(sum(P(fd1,:),2)==0))

P = Read_dmatcsr('../output/sum0_zero_row_256x256_P1.dmatcsr');
fd1 = find(~(abs(sum(P,2)-1)<0.0000000001))
fd2 = find(~(sum(P(fd1,:),2)==0))

%=========================================================
P = Read_dmatcsr('../output/fdm2d_9P_bnd_16x16_P1.dmatcsr');
s16 = sum(P,2);

P = Read_dmatcsr('../output/fdm2d_9P_bnd_64x64_P1.dmatcsr');
s64 = sum(P,2);

P = Read_dmatcsr('../output/fdm2d_9P_bnd_256x256_P1.dmatcsr');
s_256 = sum(P,2);

%=========================================================
P = Read_dmatcsr('../output/strange_16x16_P1.dmatcsr');
ss16 = sum(P,2);

