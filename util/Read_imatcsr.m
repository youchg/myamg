function A = Read_imatcsr(filename)
fid = fopen(filename);
nr = fscanf(fid, '%d', 1);
nc = fscanf(fid, '%d', 1);
nn = fscanf(fid, '%d', 1);
ia  = fscanf(fid, '%d', nr+1)';
ja  = fscanf(fid, '%d', nn)';
va  = fscanf(fid, '%d', nn)';
A   = zeros(nr, nc);
k=1;
for i=1:nr
    for j=(ia(i)+1):ia(i+1)
        A(i,ja(j)+1) = va(k);
        k=k+1;
    end
end
fclose(fid);
