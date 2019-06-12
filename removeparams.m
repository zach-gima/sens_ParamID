function [ Arem ] = removeparams( A, par)
%removes the rows and colum correspondings to the indicies in "par"
Arem=A;
for i=abs(sort(-(abs(par))))
Arem=removeparam(Arem,i);
end
end


