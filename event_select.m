%% event_select: computes the normalized sensitivity (S^T*S) content for multpiple discrete set of data (called events)
% then determines which subset of events has the highest sensitivity
% content (subset size determined by user)

function [opt_event_idx] = event_select(ID_p,STS_norm)
    
    C = STS_norm;

    %% Run CVX routine to select events with highest sensitivity
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
    
    cvx_event_idx = round(ID_p.event_budget*x); % generates binary vector indicating which events are selected; typically not an integer so we round
    opt_event_idx = find(cvx_event_idx);
    
%     Num_exp = ceil(m*x);
    obj_mat = 0;
    for k= 1: ID_p.num_events
        obj_mat = obj_mat + cvx_event_idx(k)*C{k};
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

end