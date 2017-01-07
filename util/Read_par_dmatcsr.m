function A = Read_par_dmatcsr(filename, nrank)

Apart = cell(nrank, 1);
for i = 1 : nrank
    Apart{i} = Read_dmatcsr([filename, '.par.dat.rank', num2str(i-1, '%02d')]);
end

A = [];
for i = 1 : nrank
    A = [A; Apart{i}];
end
