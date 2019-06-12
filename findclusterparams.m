%%% Adapted from Bosch Year 2 Project (function originally by Dylan Kato)

function [ removed ] = findclusterparams(STS_norm_diag,corr_coeff_matrix,ID_p)
    %% Parse Inputs
    Np = ID_p.np;
    threshold = ID_p.collinearity_thresh;
    
    %% Implement clustering logic and parameter elimination
    %This function takes the values in STSnorm and a given threshold and it
    %will find what parameters must be removed
    minSTSnorm = elementwisecellmin(corr_coeff_matrix,Np)-diag(ones(Np,1));
    removed=[];
    
    % iterate through collinearity values until the max collinearity value
    % in the sensitivity matrix is below the threshold
    while true
        % find the maximum minimum sensitivity element excluding those removed
        maxS = maxexc(minSTSnorm,removed);
        
        % check if this element is above the threshold
        % if not, end
        if maxS<threshold
            break
        end
        % if so, procede

        % find the indexes of this element
        [row,col] = find(minSTSnorm==maxS,1);
        % compare their sensitivities
        %%% ZTG fix this so that it works when only using one experiment.
        %%% look at avgSTS and Sens_mag
        senscomp = [STS_norm_diag(row),STS_norm_diag(col)];
        % save which is higher and which is lower
        if senscomp(1) >= senscomp(2)
            rem=col;
        end
        if senscomp(1) < senscomp(2)
            rem=row;
        end
        % append the lower to a list of removed elements
        removed = [removed,rem];
        % repeat for appended matrix
    end
    
    removed = sort(removed);
end

