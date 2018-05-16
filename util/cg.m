function x = cg(A, b, x0, tol, max_iter, sol)
format long
x = x0;
r = b - A * x;
rn = norm(r);
i = 0;
tmp = 9999;
prt = [];
while rn>tol && i<max_iter
    i = i+1;
    if(i == 1)
        status = 0;
        p = r;
    else
        beta = rn*rn / (tmp*tmp);
        p = p*beta + r;
    end
    
    Ap = A*p;
    alpha = rn*rn / (p'*Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    tmp = rn;
    rn = norm(r);
    error = norm(sol-x);
    prt = [prt; i, rn, error, rn/tmp];
end
fprintf(1, '%3d, %20.15f, %20.15f, %18.15f\n', prt');
