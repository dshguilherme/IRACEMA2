function KnotSpanLinear = FindSpanLinear(n, p, u, U)

%implementa��o para determinar o �ndice do knot span por uma busca linear
%entrada: n, p, u, U 
%U(knot vector), u valor arbitr�rio, p ordem e n = m-p-1
%sa�da �ndice do knot span

%m = length(U);

    KnotSpanLinear = n;
    if u == U(n+2)
       return
    end
    
    for j=1:n+1
        if(u >= U(j) && u < U(j+1))
            KnotSpanLinear = j -1;
            return
        end
    end
end

    