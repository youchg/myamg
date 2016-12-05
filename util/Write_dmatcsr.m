function status = Write_dmatcsr(A, filename)
status = 0;
fid = fopen(filename, 'w');

[nr, nc] = size(A);
fprintf(fid, '%d\n', nr);
fprintf(fid, '%d\n', nc);

nn = 0;
ia = zeros(1, nr+1);
for i = 1:nr
    for j = 1:nc
        if A(i, j) ~= 0
            nn = nn+1;
            ia(i+1) = ia(i+1) + 1;
        end
    end
end
fprintf(fid, '%d\n\n', nn);

fprintf(fid, '%d\n', ia(1));
for i = 2:nr+1
    ia(i) = ia(i) + ia(i-1);
    fprintf(fid, '%d\n', ia(i));
end
fprintf(fid, '\n');

for i = 1:nr
    for j = 1:nc
        if A(i, j) ~= 0
            fprintf(fid, '%d\n', j-1);
        end
    end
end
fprintf(fid, '\n');

for i = 1:nr
    for j = 1:nc
        if A(i, j) ~= 0
            fprintf(fid, '%15.12f\n', A(i, j));
        end
    end
end

fclose(fid);
status = 1;
