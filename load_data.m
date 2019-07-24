function data = load_data(p,input_folder)

    %%% Set Input Profile 
%     load(input_filename,'Time_exp','V_LM_CELL_sim','Current_exp','T_amb_sim')
%     Time_exp_cell = Time_exp;
%     Current_exp_cell = Current_exp;
%     clear Time_exp Current_exp
%       
%     filename = cell(6,1);
%     filename{1} = 'DC1_batt_ObsData.mat';
%     filename{2} = 'DC2_batt_ObsData.mat';
%     filename{3} = 'LA92x2_batt_ObsData.mat';
%     filename{4} = 'SC04x4_batt_ObsData.mat';
%     filename{5} = 'UDDSx2_batt_ObsData.mat';
%     filename{6} = 'US06x3_batt_ObsData.mat';
%     for ii = 1:6
%         Time = Time_exp_cell{ii};
%         Current = Current_exp_cell{ii};
%         Voltage = V_LM_CELL_sim{ii};
%         T_amb = T_amb_sim{ii};
%         save(filename{ii},'Time','Current','Voltage','T_amb');
%     end
    
    input_files = dir(input_folder);
    % on Mac, use line below to ignore '.','..', and '.DS_Store' files that
    % are loaded into r 
    input_files=input_files(~ismember({input_files.name},{'.','..','.DS_Store'}));
    num_inputs = length(input_files);
    
    data(num_inputs) = struct();
    for ii = 1:num_inputs
        input_filename = input_files(ii).name;

        load(strcat(input_folder,input_filename),'Time','Current','Voltage','T_amb','States')

        % See if Rc value loaded/specified for the input profile being loaded;
        % otherwise do nothing
        try
            p.R_c = Rc;
        catch
            fprintf('No Rc value specified for %s profile. \n',input_filename)
        end

        % Data structure with time,current, initial condition
        data(ii).cycle_name = input_filename;
        data(ii).time = Time;
%         data(ii).V_exp = Voltage;
%         data(ii).V0 = Voltage(1);
        data(ii).V0 = 3.768907470798727; % SOC = 60%; ZTG CHANGE 2019-7-22 COME BACK AND ALTER IF CHANGING INITIAL SOC
        data(ii).T_amb = T_amb;
%         data(ii).states_true = States;

        % % % Current | Positive <=> Discharge, Negative <=> Charge
        % % % Opposite convention outside the models
        data(ii).cur = Current; % not normalized by electrode area for SPMeT (but is for DFN)
        
        clear Time Voltage T_amb Current States
    end
end