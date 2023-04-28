function B = elementStiffness(dR,d)
B = zeros((d-1)*3,length(dR)*d);
switch d
    case 3
        for i=1:length(dR)
            B(:,1+(i-1)*d:3+(i-1)*d) = [dR(i,1) 0 0;
                0 dR(i,2) 0;
                0 0 dR(i,3);
                0 dR(i,3) dR(i,2);
                dR(i,3) 0 dR(i,1);
                dR(i,2) dR(i,1) 0];
        end
    case 2
        for i=1:length(dR)
            B(:,1+(i-1)*d:2+(i-1)*d) = [dR(i,1) 0;
                0 dR(i,2);
                dR(i,2) dR(i,1)];
        end

end