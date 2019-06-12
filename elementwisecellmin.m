function [ minmat ] = elementwisecellmin( STSnorm,Np )
    %takes the elementwise minimum of the cell array STS
    minn = ones(Np,Np,2);
    Z = STSnorm;
    minn(:,:,2)=Z;
    minn(:,:,1)=min(minn,[],3);
    minmat = minn(:,:,1);
end

