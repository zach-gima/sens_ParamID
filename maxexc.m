function [ maxx ] = maxexc( A,b )
    %finds the max of A excluding rows and columns in b
    b=sort(-b);

    for i=abs(b)
        A=removeparam(A,i);
    end
    maxx=max(max(A));
end

