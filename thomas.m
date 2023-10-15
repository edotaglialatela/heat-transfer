% Funzione per risolvere un sistema tridiagonale con il metodo di Thomas
%[A.Quarteroni - Matematica Numerica]
function x = thomas(A, b)
    n = length(b);
    x = zeros(n, 1);
    c = diag(A);
    d = diag(A, 1);
    e = diag(A, -1);
    
    % Fattorizzazione LU
    for k = 2:n
        factor = e(k-1) / c(k-1);
        c(k) = c(k) - factor * d(k-1);
        b(k) = b(k) - factor * b(k-1);
    end
    
    % Sostituzione in avanti
    x(n) = b(n) / c(n);
    for k = n-1:-1:1
        x(k) = (b(k) - d(k) * x(k+1)) / c(k);
    end
end