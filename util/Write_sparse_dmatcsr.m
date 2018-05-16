function status = Write_sparse_dmatcsr(A, filename)
status = 0;
fid    = fopen(filename, 'w');

%====================================================================
[nr, nc] = size(A);
fprintf(fid, '%d\n', nr);
fprintf(fid, '%d\n', nc);
fprintf(1,   '\nWrite matrix: (%d x %d)\n\n', nr, nc);

%====================================================================
nn = nnz(A);
fprintf(1,   'The number of nonzeros is %d\n\n', nn);
fprintf(fid, '%d\n\n', nn);

%====================================================================
percent_divide = nr / 20;
percent        = 1;

%====================================================================
S  = spones(A);
ia = [0; full(sum(S, 2))];

percent = 1;
fprintf(1,   'Write matrix->ia...\n');
fprintf(1,   '%d%% ', 0);
fprintf(fid, '%d\n', ia(1));

nrp1 = nr + 1;
for i = 2:nrp1
    if i > percent * percent_divide
	fprintf(1, '%d%% ', percent*5);
	percent = percent + 1;
    end

    ia(i) = ia(i) + ia(i-1);
    fprintf(fid, '%d\n', ia(i));
end
fprintf(1,   '\n');
fprintf(fid, '\n');

%====================================================================
percent = 1;
fprintf(1, 'Write matrix->ja...\n');
fprintf(1, '%d%% ', 0);

for i = 1:nr
    if i > percent * percent_divide
	fprintf(1, '%d%% ', percent*5);
	percent = percent + 1;
    end

    [~, ja_row] = find(A(i, :));
    len_ja_row = length(ja_row);
    for j = 1:len_ja_row
	fprintf(fid, '%d\n', ja_row(j)-1);
    end
end
fprintf(1,   '\n');
fprintf(fid, '\n');

%====================================================================
va     = nonzeros(A');
len_va = length(va);

percent_divide = len_va / 20;
percent        = 1;

fprintf(1, 'Write matrix->va...\n');
fprintf(1, '%d%% ', 0);
for i = 1:len_va
    if i > percent * percent_divide
	fprintf(1, '%d%% ', percent*5);
	percent = percent + 1;
    end

    fprintf(fid, '%15.12f\n', va(i));
end

%====================================================================
fclose(fid);
status = 1;
