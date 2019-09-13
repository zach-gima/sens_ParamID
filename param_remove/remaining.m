function [ remainins ] = remaining( Np, removed)
%removes the rows and colum correspondings to the indicies in "par"
remainins = [1:Np];
remainins=removecols(remainins,removed);
end