function [g, dg] = topologicalInterp(c, inter_geo)
    if c < 0.1
        g = 0.1;
        dg = 0;
    elseif c > 0.945
        g = 1;
        dg = 0;
    else
        g = 0.1 +(1-0.1)*inter_geo.eval_point(c);
        g = g(1);
        dg = 0.9*inter_geo.parameter_derivative([],c);
        dg = dg(1);
    end
    g = g^3;
    dg = 3*(g^2)*dg;
end