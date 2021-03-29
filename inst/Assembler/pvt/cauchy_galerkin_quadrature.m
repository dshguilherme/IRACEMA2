function [qpoints, qweights] = cauchy_galerkin_quadrature(n)
switch n
    case 3
        tmp = 1/sqrt(3);
        qpoints(1) = -tmp;
        qpoints(2) = tmp;
    case 4
        qpoints(1) = -1;
        qpoints(2) = 0;
        qpoints(3) = 1;
    case 5
        tmp = sqrt(225 -30*sqrt(30))/15;
        qpoints(1) = -tmp;
        qpoints(2) = tmp;
    case 6
        qpoints(1) = -1;
        qpoints(2) = 0;
        qpoints(3) = 1;
    case 7
        tmp = 0.504918567512;
        qpoints(1) = -tmp;
        qpoints(2) = tmp;
    case 8
        qpoints(1) = -1;
        qpoints(2) = 0;
        qpoints(3) = 1;
    otherwise
end
qweights = ones(size(qpoints));
        
end