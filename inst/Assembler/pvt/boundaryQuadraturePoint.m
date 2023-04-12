function q = boundaryQuadraturePoint(qp, side)
    if mod(side,2) == 0
        b = 1;
    else
        b = -1;
    end
    
    idx = round(side/2);
    [sz1, sz2] = size(qp);
    n_idx = setdiff(sz2, idx);
    
    q = zeros(, sz2);
    
    qp(:,idx) = b;
    qp(:,n_idx) = points;

end