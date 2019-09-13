function [ Arem ] = removecol( A,col )
%removes a column of a matrix
Arem = [A(:,1:col-1),A(:,col+1:end)];
end

