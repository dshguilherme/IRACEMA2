function f = cantileverForceAdjust(x, L, h1, h2, intensity)
fx = 0;
if x(1) > L && (x(2) < h1 || x(2) > h2 )
    fy = -intensity;
else
    fy = 0;
end
f = [fx fy 0];
end


%----------------------------------------           
%                                       |   h1       
%                                       |___|_____      
%                                       |        | x(2) < h1 || x(2) > h2  
%                                       |________|
%                                       |   |    | 
%                                       |   h2   |
%----------------------------------------   |    |