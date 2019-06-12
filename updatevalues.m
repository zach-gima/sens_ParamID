function [ A, Sens_orth ] = updatevalues( sens, removed, NT, Np, r)
    %this will update the important values after removing the parameters in
    %"removed"
    
    Np_rem=Np-length(removed);
    A=zeros(Np_rem,length(sens));
    Sens_orth=zeros(Np_rem,length(sens));
    mu = 10^-4;
    for i=1:length(sens)
        name=r(i).name; % removed + 2 and updated r abovZ
        index=str2double(name(1:end-4));
        S3 = sens{i};
        S3 = S3(:,remaining(Np,removed)); %Note this should be the original Np size before removing)

        [~,R,E] = qr(S3,0);
        % Extract Diagonal
        Di = diag(R);
        for j = 1:Np_rem
             D(j) = Di(find(E == j));
        end

        % apply sensitivity threshold; if parameters not above noise
        % sensitivitiy threshold, then A will store a 0; otherwise the
        % experiment index
        A(:,i) = index*(abs(D)>(mu*sqrt(NT(i))));
        Sens_orth(:,i)=abs(D);
    end
end

