%% cvx_ParamID function
% ZTG Note 2019-6-6: previous implementation of batch select, although some
% variables have already been changed so this is a partially altered
% function
function [es_out] = batch_select(ID_p,data)
    %% Parse Inputs    
    t = data.time;
    I = data.cur;
    V = data.V_exp;
%     sens = data.sens;
    sens_norm = data.sens_norm;
    
    %% Split cycle data into batches of pre-determined size and store in cell
%     num_batches = ceil(length(I)/batch_size); 
    
    % This generates a length(I) x 1 vector, where every batch_size segment of data has the corresponding batch number. 
    % e.g. For batch size = 50, batch_idx(1:50) = 1; batch_idx(51:100) = 2;
    event_idx = ceil((1:length(I))/batch_size); 
    
    C = cell(ID_p.num_events,1);
    t_event = cell(ID_p.num_events,1);
    I_event = cell(ID_p.num_events,1);
    V_event = cell(ID_p.num_events,1);
    V0_event = cell(ID_p.num_events,1);
    STS_event = cell(ID_p.num_events,1);
    STSmag_event = cell(ID_p.num_events,1);
    T_amb_event = cell(ID_p.num_events,1);
    STSnorm_batch = cell(ID_p.num_events,1);
%     c_s_n_batch = cell(ID_p.num_events,1);
    
    % Set covariance
%     Q = 1; % simulation case
    
    % Batch it up baby!!
    for ii = 1:ID_p.num_events
        % Input Data
        t_event{ii} = t(event_idx == ii);
        I_event{ii} = I(event_idx == ii);
%         dV_batch = dV(batch_idx == ii);
        V_event{ii} = V(event_idx == ii);
        V0_event{ii} = V_event{ii}(1);
%         c_s_n_batch{ii} = c_nx(batch_idx == ii,:);
        T_amb_event{ii} = data.T_amb;

        % Compute Sensitivity Content
        sens_batch = sens_norm(event_idx == ii,:);
        STS_event{ii} = sens_batch'*sens_batch;
        STSmag_event{ii} = diag(STS_event{ii});
        ttt=inv(diag(sqrt(diag(STS_event{ii}))));
        STSnorm_batch{ii}=abs(ttt*STS_event{ii}*ttt);  
        
        C{ii} = STS_event{ii};
    end

    %% Run CVX routine to select batches with highest sensitivity
    disp('Starting cvx batch select');

    cvx_begin
        cvx_solver SeDuMi

        variable x(ID_p.num_events)

        bigOne = ones(1,ID_p.num_events);
        expression M

        for k = 1: ID_p.num_events
            M = M + x(k)*C{k};
        end

        % D-optimality condition
        maximize(det_rootn(M))

        subject to
%             m*cost'*x <= B
            bigOne*x == 1; %
            0 <= x <= 1/ID_p.event_budget;
    cvx_end
    
    cvx_batch_idx = round(ID_p.event_budget*x);
    opt_batch_num = find(cvx_batch_idx);
    
%     Num_exp = ceil(m*x);
    obj_mat = 0;
    for k= 1: ID_p.num_events
        obj_mat = obj_mat + cvx_batch_idx(k)*C{k};
    end
    fprintf('Determinant of STS: %e \n',det(obj_mat));

%     % Plot stuff
%     figure()
%     bar(cvx_batch_idx)
%     xlabel('Batch ID #')
%     xlim([1,ID_p.num_events])
%     ylabel('# of Times Chosen')
%     set(gca,'FontSize',16)
%     
%     % Superimpose data selected on top of original current profile
%     % figure('Position', [100 100 900 700])
%     figure('Position', [100 100 900 1600])
%     plot(t,I,'Color',[1, 0.5, 0],'LineWidth',3); %plot full I profile
%     hold on
%     for jj = 1:length(opt_batch_num)
%         plot(t_batch{opt_batch_num(jj)},I_batch{opt_batch_num(jj)},'*','Color',[0, 0.6, 0],'MarkerSize',10,'MarkerEdgeColor','k');
%     end
%     hold off
%     title('Input Current Profile');
%     xlabel('Time (s)','FontSize',fs);
%     ylabel('Current (A)','FontSize',fs)
%     set(gca,'FontSize',fs)
%     box on
%     grid on    
%     
%         
    %% Store Outputs
    es_out.opt_batch_num = opt_batch_num;
    es_out.I_batch = I_event;
    es_out.t_batch = t_event;
    es_out.V_batch = V_event;
    es_out.V0_batch = V0_event;
    es_out.STS_batch = STS_event;
    es_out.STSmag_batch = STSmag_event;
    es_out.STSnorm_batch = STSnorm_batch;
    es_out.T_amb_batch = T_amb_event;
%     es_out.c_s_n_batch = c_s_n_batch;
end