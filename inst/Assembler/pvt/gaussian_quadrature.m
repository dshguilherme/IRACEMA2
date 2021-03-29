function [qpoints, qweights] = gaussian_quadrature(number_of_points)
    coefs = LegendrePoly(number_of_points);
    qpoints = sort(roots(coefs));
    qweights = zeros(size(qpoints));
    der = polyder(coefs);

    for i=1:numel(qpoints)
        tmp1 = polyval(der,qpoints(i));
        tmp2 = (1- qpoints(i)^2)*(tmp1^2);

        qweights(i) = 2/tmp2;
    end
end