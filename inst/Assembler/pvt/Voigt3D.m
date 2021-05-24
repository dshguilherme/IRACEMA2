function alpha = Voigt3D(i,j)
if i == j
    alpha = i;
else
    switch i
        case 1
            if j == 2
                alpha = 6;
            elseif j == 3
                alpha = 5;
            end
        case 2
            alpha = 4;
    end
end
                
end