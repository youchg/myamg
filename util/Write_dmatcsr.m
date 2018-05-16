function status = Write_dmatcsr(A, filename)
status = 0;
fid = fopen(filename, 'w');

[nr, nc] = size(A);
fprintf(fid, '%d\n', nr);
fprintf(fid, '%d\n', nc);

fprintf(1, '\nWrite matrix: (%d x %d)\n\n', nr, nc);

%====================================================================
fprintf(1, 'Find number of nonzeros...\n');
percent_divide = floor(nr / 20);
percent_divide = nr / 20;
percent = 1;
fprintf(1, '%d%% ', 0);

nn = 0;
ia = zeros(1, nr+1);
for i = 1:nr
    if i > percent * percent_divide
	fprintf(1, '%d%% ', percent*5);
	percent = percent + 1;
    end

    for j = 1:nc
        if A(i, j) ~= 0
            nn = nn+1;
            ia(i+1) = ia(i+1) + 1;
        end
    end
end
fprintf(1, '\n');

fprintf(fid, '%d\n\n', nn);
fprintf(1, 'The number of nonzeros is %d\n\n', nn);

%====================================================================
fprintf(1, 'Write matrix->ia...\n');
percent = 1;
fprintf(1, '%d%% ', 0);

fprintf(fid, '%d\n', ia(1));
for i = 2:nr+1
    if i > percent * percent_divide
	fprintf(1, '%d%% ', percent*5);
	percent = percent + 1;
    end

    ia(i) = ia(i) + ia(i-1);
    fprintf(fid, '%d\n', ia(i));
end
fprintf(1, '\n');
fprintf(fid, '\n');

%====================================================================
fprintf(1, 'Write matrix->ja...\n');
percent = 1;
fprintf(1, '%d%% ', 0);

for i = 1:nr
    if i > percent * percent_divide
	fprintf(1, '%d%% ', percent*5);
	percent = percent + 1;
    end
    for j = 1:nc
        if A(i, j) ~= 0
            fprintf(fid, '%d\n', j-1);
        end
    end
end
fprintf(1, '\n');
fprintf(fid, '\n');

%====================================================================
fprintf(1, 'Write matrix->va...\n');
percent = 1;
fprintf(1, '%d%% ', 0);

for i = 1:nr
    if i > percent * percent_divide
	fprintf(1, '%d%% ', percent*5);
	percent = percent + 1;
    end

    for j = 1:nc
        if A(i, j) ~= 0
            fprintf(fid, '%15.12f\n', A(i, j));
        end
    end
end
fprintf(1, '\n');

%====================================================================
fclose(fid);
status = 1;
