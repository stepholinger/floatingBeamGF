function c = cubic_basis(x,n,derivative)

if n == 1
    if derivative == 0
        c = (3-2*x)*x^2;
    elseif derivative == 1
        c = 6*x-6*x^2;
    elseif derivative == 2
        c = 6-12*x;
    end
elseif n == 2
    if derivative == 0
        c = -1*(1-x)*x^2;
    elseif derivative == 1
        c = -2*x+3*x^2;
    elseif derivative == 2
        c = 6*x-2;
    end
elseif n == 3
    if derivative == 0    
        c = x*(x-1)^2;
    elseif derivative == 1
        c = 3*(x^2)-4*x+1;
    elseif derivative == 2
        c = 6*x-4;
    end
elseif n == 4
    if derivative == 0
        c = 2*(x^3)-3*(x^2)+1;
    elseif derivative == 1    
        c = 6*(x^2)-6*x;
    elseif derivative == 2
        c = 12*x-6;
    end
end

end