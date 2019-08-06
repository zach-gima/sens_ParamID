function [metrics] = compute_metrics(num_batches,partial_path)
    run param/params_bounds

    %%% Pre-allocate states of interest
    % Surface concentration
    cssn_initial = [];
    cssn_final = [];
    cssp_initial = [];
    cssp_final = [];
    % Side rxn overpotential
    etas_initial = [];
    etas_final = [];    

    %Electrode overpotential
    %     etan_initial = [];
    %     etan_final = [];
    %     etap_initial  = [];
    %     etap_final  = [];

    % elyte conc. at current collectors
    ce0n_initial = [];
    ce0n_final = [];
    ce0p_initial = [];
    ce0p_final = [];

    cssn_true = [];
    cssp_true = [];
    etas_true =  [];
    ce0n_true =  [];
    ce0p_true =  [];
    
    rmse_vec = [];
    cost_evolution = cell(num_batches,1);
    
    % fmincon stuff
    fmincon_iter = zeros(num_batches,1);
    
    for batch_idx = 1:num_batches
        % Load data
        full_path = strcat(partial_path,num2str(batch_idx),'.mat');
        load(full_path);

        %% Parse Data
        opt_event_idx = ID_out.opt_event_idx;
        num_opt_events = length(opt_event_idx);
        data = ID_out.data;
        opt_data = data(opt_event_idx);

        V_sim_initial = cell2mat(ID_out.V_sim_initial(opt_event_idx));
        V_sim_final = cell2mat(ID_out.V_sim_final(opt_event_idx));
        
        states_initial = ID_out.states_initial(opt_event_idx);
        states_final = ID_out.states_final(opt_event_idx);

        theta = ID_out.theta;
        theta_iter_vec = 0:num_batches;

        fmincon_iter(batch_idx) = ID_out.fmincon_out.iterations + 1; % +1 to includ e initial condition

        %% Parse Truth & Time Data
        norm_truth_param = origin_to_norm(theta.truth,p_bounds);
        
        t_end = 0; % for concatenating the time vector data

        for zz = 1:num_opt_events
            t_cell{zz,1} = opt_data(zz).time + t_end;
            V_true_cell{zz,1} = opt_data(zz).V_exp;
            states_true{zz,1} = opt_data(zz).states_true;

            cssn_true = vertcat(cssn_true,states_true{zz,1}.cssn_sim);
            cssp_true= vertcat(cssp_true,states_true{zz,1}.cssp_sim);
            etas_true = vertcat(etas_true,states_true{zz,1}.etas_sim);
            ce0n_true = vertcat(ce0n_true,states_true{zz,1}.ce0n_sim);
            ce0p_true = vertcat(ce0p_true,states_true{zz,1}.ce0p_sim);

            t_end = t_end + opt_data(zz).time(end);
        end

        t = cell2mat(t_cell);
        V_true = cell2mat(V_true_cell);

        %% Voltage Stuff
        if batch_idx == 1
            cost_evolution{batch_idx} = ID_out.fmincon_history.fval/num_opt_events;
            rmse_vec = vertcat(rmse_vec,[cost_evolution{batch_idx}(1);cost_evolution{batch_idx}(end)]);
%             rmse_vec_2 = vertcat(rmse_vec_2,[rmse(V_sim_initial,V_true);rmse(V_sim_final,V_true)]);
        else
            cost_evolution{batch_idx} = ID_out.fmincon_history.fval/num_opt_events;
            rmse_vec = vertcat(rmse_vec,cost_evolution{batch_idx}(end));
%             rmse_vec_2 = vertcat(rmse_vec_2,rmse(V_sim_final,V_true));
        end
        
        %below if structure to handle weird way fmincon outputs were saved;
        % should be able to remove in future
%         if batch_idx == 1
%             cost_evolution{1} = ID_out.fmincon_history.fval;
%             rmse_vec = vertcat(rmse_vec,[cost_evolution{1}(1);cost_evolution{1}(end)]);
%         else
%             start_idx = fmincon_iter(batch_idx-1)+1;
%             cost_evolution{batch_idx} = ID_out.fmincon_history.fval(start_idx:end);
%             rmse_vec = vertcat(rmse_vec,cost_evolution{batch_idx}(end));
%         end

        %% Internal State RMSE

        % For first batch, compute initial rmse before ID
        if batch_idx == 1
            for jj = 1:num_opt_events
                cssn_initial = vertcat(cssn_initial,states_initial{jj}.cssn_sim);
                cssp_initial = vertcat(cssp_initial,states_initial{jj}.cssp_sim);
                ce0p_initial = vertcat(ce0p_initial,states_initial{jj}.ce0p_sim);
                ce0n_initial = vertcat(ce0n_initial,states_initial{jj}.ce0n_sim);
                etas_initial = vertcat(etas_initial,states_initial{jj}.etas_sim);
        %         etan_initial = vertcat(etan_initial,states_initial{jj}.etan_sim);
        %         etap_initial = vertcat(etap_initial,states_initial{jj}.etap_sim);
            end

            cssn_rmse(1) = rmse(cssn_true,cssn_initial);
            cssp_rmse(1) = rmse(cssp_true,cssp_initial);
            etas_rmse(1) = rmse(etas_true,etas_initial);
            ce0n_rmse(1) = rmse(ce0n_true,ce0n_initial);
            ce0p_rmse(1) = rmse(ce0p_true,ce0p_initial);
    %       etan_rmse(1) = rmse(etan_true,etan_initial);
    %       etap_rmse(1) = rmse(etap_true,etap_initial);
            norm_initial_param = origin_to_norm(theta.guess,p_bounds);
            norm_param_dist(1) = norm(norm_truth_param-norm_initial_param);
            per_param_norm_error{1} = (norm_truth_param-norm_initial_param)./norm_truth_param*100;
        end

        for jj = 1:num_opt_events
            cssn_final = vertcat(cssn_final,states_final{jj}.cssn_sim);
            cssp_final = vertcat(cssp_final,states_final{jj}.cssp_sim);        
            ce0n_final = vertcat(ce0n_final,states_final{jj}.ce0n_sim);
            ce0p_final = vertcat(ce0p_final,states_final{jj}.ce0p_sim);
            etas_final = vertcat(etas_final,states_final{jj}.etas_sim);
    %         etan_final = vertcat(etan_final,states_final{jj}.etan_sim);
    %         etap_final = vertcat(etap_final,states_final{jj}.etap_sim);
        end

        %%% css RMSE and fits
        cssn_rmse(batch_idx+1) = rmse(cssn_true,cssn_final);
        cssp_rmse(batch_idx+1) = rmse(cssp_true,cssp_final);

        %%% electrode overpotential RMSE and fits    
    %     etan_rmse(batch_idx+1) = rmse(etan_true,etan_final);
    %     etap_rmse(batch_idx+1) = rmse(etap_true,etap_final);

        %%% sid rxn overpotential RMSE and fits
        etas_rmse(batch_idx+1) = rmse(etas_true,etas_final);

        %%% electrode concentration @ current collectors
        ce0n_rmse(batch_idx+1) = rmse(ce0n_true,ce0n_final);
        ce0p_rmse(batch_idx+1) = rmse(ce0p_true,ce0p_final);

        %% normalized parameter distance
        norm_ID_param = origin_to_norm(theta.final_ID,p_bounds);
        norm_param_dist(batch_idx+1) = norm(norm_truth_param-norm_ID_param);
        per_param_norm_error{batch_idx+1} = (norm_truth_param-norm_ID_param)./norm_truth_param*100;
        
        clear ID_out
    end
    
    % save and output data
    metrics.cssn_rmse = cssn_rmse;
    metrics.cssp_rmse = cssp_rmse;
    metrics.etas_rmse = etas_rmse;
    metrics.ce0n_rmse = ce0n_rmse;
    metrics.ce0p_rmse = ce0p_rmse;
    metrics.norm_param_dist = norm_param_dist;
    metrics.fmincon_iter = fmincon_iter;
%     metrics.cost_evolution = cost_evolution;
    metrics.rmse_vec = rmse_vec;
    metrics.theta_iter_vec = theta_iter_vec;
    metrics.per_param_norm_error = per_param_norm_error;
    
end