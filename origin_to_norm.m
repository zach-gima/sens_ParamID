%% Function for converting parameters from nominal values to normalized values
% By: Zach Gima 2019-6-28

% want to only return the parameters being identified

function [norm_param] = origin_to_norm(origin_param, bounds) %, paramID_idx)
    % Note that this is currently written for a parameter set where each
    % parameter as the same normalization scheme; some parameters do not
    % use this scheme so will have to rewrite / expand functionality
%     norm_param = (origin_param(paramID_idx) - bounds.min(paramID_idx)) ./ (bounds.max(paramID_idx) - bounds.min(paramID_idx));
    norm_param = (origin_param - bounds.min) ./ (bounds.max - bounds.min);
end

