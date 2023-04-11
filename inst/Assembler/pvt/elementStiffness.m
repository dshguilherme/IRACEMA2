function B = elementStiffness(dR, d, nel_dof)
    B = zeros(d, nel_dof*d);
    for i=1:d
        B(i,i:d:end) = dR(:,i);
    end
    if d == 2
        B(3,1:d:end) = dR(:,2);
        B(3,2:d:end) = dR(:,1);
    elseif d == 3
        B(4,2:d:end) = dR(:,3);
        B(4,3:d:end) = dR(:,2);
        B(5,1:d:end) = dR(:,3);
        B(5,3:d:end) = dR(:,1);
        B(6,1:d:end) = dR(:,2);
        B(6,2:d:end) = dR(:,1);
    end
end