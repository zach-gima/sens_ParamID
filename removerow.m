function [ Arem ] = removerow( A,row )
%removes a row of a matrix
AremT=removecol(A',row);
Arem=AremT';
end