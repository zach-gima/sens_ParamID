function [ Arem ] = removeparam( A, par)
%removes the row and colum corresponding to the index "par"
Arem=removecol(A,par);
Arem=removerow(Arem,par);
end

