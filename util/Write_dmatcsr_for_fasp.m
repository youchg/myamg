function status = Write_dmatcsr_for_fasp(A, filename)
status = 0;
fid = fopen(filename, 'w');

[nr, nc] = size(A);
fprintf(fid, '%d\n', nr);

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

%progress = 1

fprintf(fid, '%d\n', ia(1)+1);
for i = 2:nr+1
    ia(i) = ia(i) + ia(i-1);
    fprintf(fid, '%d\n', ia(i)+1);
end
fprintf(fid, '\n');

%progress = 2

for i = 1:nr
    for j = 1:nc
        if A(i, j) ~= 0
            fprintf(fid, '%d\n', j);
        end
    end
end
fprintf(fid, '\n');

for i = 1:nr
    for j = 1:nc
        if A(i, j) ~= 0
            fprintf(fid, '%e\n', A(i, j));
        end
    end
end

%progress = 3

fclose(fid);
status = 1;
