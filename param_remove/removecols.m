function [ Arem ] = removecols( A, par)
%removes the rows and colum correspondings to the indicies in "par"
Arem=A;
for i=abs(sort(-(abs(par))))
Arem=removecol(Arem,i);
end
end