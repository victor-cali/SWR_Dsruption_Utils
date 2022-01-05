classdef RipplesAnalyzerMatFile < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        Panel                         matlab.ui.container.Panel
        EventDropDown                 matlab.ui.control.DropDown
        EventDropDownLabel            matlab.ui.control.Label
        BoutDropDown                  matlab.ui.control.DropDown
        BoutDropDownLabel             matlab.ui.control.Label
        TotalHPCripplesLabel          matlab.ui.control.Label
        TotalPFCripplesLabel          matlab.ui.control.Label
        BoutMoreRipplesLabel          matlab.ui.control.Label
        EpochDropDown                 matlab.ui.control.DropDown
        EpochDropDownLabel            matlab.ui.control.Label
        WindowLengthEditField         matlab.ui.control.NumericEditField
        WindowlengthmsEditFieldLabel  matlab.ui.control.Label
        BrainAreaButtonGroup          matlab.ui.container.ButtonGroup
        HPCButton                     matlab.ui.control.ToggleButton
        PFCButton                     matlab.ui.control.ToggleButton
        plot_UIAxes                   matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        Epochs = struct; % Description
        ss=3;
        fs = 1000;
        D1;
        D2;
        % Dataset path
        dataset_path = 'D:\Dev\MATLAB\GenzelLab\rat5';
    end
    
    methods (Access = private)
        
        function calibrate(app)
            
            power_line_noise = 50 : 50 : 300;
            
            G=getfolder();
            % Structs with all signals to be analyzed
            for j=1:length(G) % Iterate across study days.
                cd(app.dataset_path)
                cd(G{j})
                prepost=getfolder();
                for i=1: length(prepost)
                    field = strcat(G{j}(38:end),'__',prepost{i}(21:end));
                    cd(app.dataset_path)
                    cd(G{j})
                    cd(prepost{i}) 
                    %Read brain areas and load states file.
                    HPC=dir(strcat('*','HPC','*.mat'));
                    HPC=HPC.name;
                    HPC=load(HPC);
                    HPC=getfield(HPC,'HPC_ripple');
                    HPC=HPC.*(0.195);
                    HPC = notch_filter(HPC, app.fs, power_line_noise);
                    % Cortical
                    Cortex=dir(strcat('*','PFC','*.mat'));
                    Cortex=Cortex.name;
                    Cortex=load(Cortex);
                    Cortex=getfield(Cortex,'PFC');
                    Cortex=Cortex.*(0.195);
                    Cortex = notch_filter(Cortex, app.fs, power_line_noise);
                    %Load sleep scoring
                    A = dir('*states*.mat');
                    A={A.name};
                    if  ~isempty(A)
                           cellfun(@load,A);
                    else
                          error('No Scoring found')    
                    end
                    [sd_swr]=find_std(HPC,Cortex,states,app.ss);
                    Sd_Swr.sd5_hpc_co(j,i)=sd_swr.sd5_hpc_co;
                    Sd_Swr.sd5_pfc_co(j,i)=sd_swr.sd5_pfc_co;
                    
                    app.Epochs.(field) = struct('HPC',HPC,'PFC',Cortex,'States',states);
                    
                end
            end
            
            thresholds_perday_hpc=mean(Sd_Swr.sd5_hpc_co,2);
            thresholds_perday_cortex=mean(Sd_Swr.sd5_pfc_co,2);
            
            tr(1)=mean(thresholds_perday_hpc);
            tr(2)=mean(thresholds_perday_cortex);
            
            offset1= 5;
            offset2= 0;
            
            app.D1 = round(tr(1) + offset1);
            app.D2 = round(tr(2) + offset2);
            
            app.EpochDropDown.Items = fieldnames(app.Epochs);
            
        end
        
        function update_values(app)
            
            cla(app.plot_UIAxes)
            title(app.plot_UIAxes, 'Ripple Analyzer')
            xlabel(app.plot_UIAxes, 'Time (s)')
            app.plot_UIAxes.YLim = [30 360];
            app.plot_UIAxes.YTick = [100 150 220 290];
            app.plot_UIAxes.YTickLabel = {'HPC'; 'PFC'; 'HPC (Bp)'; 'PFC (Bp)'};

            % PROPERTY
            epoch = app.Epochs.(app.EpochDropDown.Value);
            
            [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc] = swr_check_thr(epoch.HPC,epoch.PFC,epoch.States,app.ss,app.D1,app.D2);
            
            % PROPERTY
            selectedButton = app.BrainAreaButtonGroup.SelectedObject;
            BR = selectedButton.Text;
            
            if strcmp(BR,'PFC')
                max_length=cellfun(@length,swr_pfc(:,1));
                N=max_length==max(max_length);
                bouts = find(max_length);
                app.BoutDropDown.Items = cellfun(@num2str,num2cell(bouts),'uni',false);
                n=str2num(app.BoutDropDown.Value);
                sn=swr_pfc{n,3};
            else
                max_length=cellfun(@length,swr_hpc(:,1));
                N=max_length==max(max_length);
                bouts = find(max_length);
                app.BoutDropDown.Items = cellfun(@num2str,num2cell(bouts),'uni',false);
                n=str2num(app.BoutDropDown.Value);
                sn=swr_hpc{n,3};
            end
            app.EventDropDown.Items = cellfun(@num2str,num2cell(1:length(sn)),'uni',false);
            
            app.BoutMoreRipplesLabel.Text = strcat('Bout with more ripples: ', num2str(length(sn)));
            app.TotalPFCripplesLabel.Text = strcat('Total PFC: ', num2str(sum(s_pfc)));
            app.TotalHPCripplesLabel.Text = strcat('Total HPC: ', num2str(sum(s_hpc)));
            
            hpc=V_hpc{N};
            pfc=V_pfc{N};
            hpc2=signal2_hpc{N};
            pfc2=signal2_pfc{N};
            
            if length(n)>1
                'Multiple epochs with same number of events'
            end
            
            % PROPERTY
            n = str2num(app.BoutDropDown.Value);
            % PROPERTY
            win_len = app.WindowLengthEditField.Value;
            % PROPERTY
            event = str2num(app.EventDropDown.Value);
            
            increment = 10;
            plot(app.plot_UIAxes,(1:length(hpc))./app.fs, increment.*zscore(hpc)+100, 'Color', 'black')
            hold(app.plot_UIAxes,'on')
            plot(app.plot_UIAxes,(1:length(pfc))./app.fs, increment.*zscore(pfc)+150, 'Color', 'black')
            plot(app.plot_UIAxes,(1:length(hpc2))./app.fs, increment.*zscore(hpc2)+220, 'Color', 'black')
            plot(app.plot_UIAxes,(1:length(pfc2))./app.fs, increment.*zscore(pfc2)+290, 'Color', 'black')
            
            if isempty(sn)
                errordlg('No Events found','Error');
            end
            
            if ~ isempty(swr_hpc{n})
                stem(app.plot_UIAxes,[swr_hpc{n,3}], ones(length([swr_hpc{n}]),1).*50, 'Color', 'blue') %(HPC)
            end
            if ~ isempty(swr_pfc{n})
                stem(app.plot_UIAxes,[swr_pfc{n,3}], ones(length([swr_pfc{n}]),1).*50, 'Color', 'red') %Seconds (Cortex)
            end
            
            app.plot_UIAxes.XLim = [(sn(event)-win_len/app.fs), (sn(event)+win_len/app.fs)];

            hold(app.plot_UIAxes,'off')
            
            msgbox('UPDATED');
            
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            
            %Load ADRITOOLS
            addpath(genpath('D:\Dev\MATLAB\GenzelLab\ADRITOOLS'))
            
            %Load CorticoHippocampal
            addpath(genpath('D:\Dev\MATLAB\GenzelLab\CorticoHippocampal'))
            
            %Load RGS14 github
            addpath(genpath('D:\Dev\MATLAB\GenzelLab\LFP_RGS14'))
            
            cd(app.dataset_path)
            
            app.calibrate()
            app.update_values()
            
        end

        % Value changed function: WindowLengthEditField
        function WindowLengthEditFieldValueChanged(app, event)
            app.update_values()
        end

        % Value changed function: EpochDropDown
        function EpochDropDownValueChanged(app, event)
            app.update_values()
        end

        % Value changed function: BoutDropDown
        function BoutDropDownValueChanged(app, event)
            app.update_values()
        end

        % Value changed function: EventDropDown
        function EventDropDownValueChanged(app, event)
            app.update_values()
        end

        % Selection changed function: BrainAreaButtonGroup
        function BrainAreaButtonGroupSelectionChanged(app, event)
            app.update_values()
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {'1x', 147};
            app.GridLayout.Padding = [21.5 10 21.5 10];

            % Create plot_UIAxes
            app.plot_UIAxes = uiaxes(app.GridLayout);
            title(app.plot_UIAxes, 'Ripple Analyzer')
            xlabel(app.plot_UIAxes, 'Time (s)')
            app.plot_UIAxes.YLim = [30 360];
            app.plot_UIAxes.YTick = [100 150 220 290];
            app.plot_UIAxes.YTickLabel = {'HPC'; 'PFC'; 'HPC (Bp)'; 'PFC (Bp)'};
            app.plot_UIAxes.Layout.Row = 1;
            app.plot_UIAxes.Layout.Column = 1;

            % Create Panel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Title = 'Panel';
            app.Panel.Layout.Row = 2;
            app.Panel.Layout.Column = 1;

            % Create BrainAreaButtonGroup
            app.BrainAreaButtonGroup = uibuttongroup(app.Panel);
            app.BrainAreaButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @BrainAreaButtonGroupSelectionChanged, true);
            app.BrainAreaButtonGroup.Title = 'Brain Area';
            app.BrainAreaButtonGroup.Position = [269 10 100 106];

            % Create PFCButton
            app.PFCButton = uitogglebutton(app.BrainAreaButtonGroup);
            app.PFCButton.Text = 'PFC';
            app.PFCButton.Position = [27 52 45 22];
            app.PFCButton.Value = true;

            % Create HPCButton
            app.HPCButton = uitogglebutton(app.BrainAreaButtonGroup);
            app.HPCButton.Text = 'HPC';
            app.HPCButton.Position = [27 11 45 22];

            % Create WindowlengthmsEditFieldLabel
            app.WindowlengthmsEditFieldLabel = uilabel(app.Panel);
            app.WindowlengthmsEditFieldLabel.HorizontalAlignment = 'right';
            app.WindowlengthmsEditFieldLabel.Position = [6 94 112 22];
            app.WindowlengthmsEditFieldLabel.Text = 'Window length (ms)';

            % Create WindowLengthEditField
            app.WindowLengthEditField = uieditfield(app.Panel, 'numeric');
            app.WindowLengthEditField.ValueChangedFcn = createCallbackFcn(app, @WindowLengthEditFieldValueChanged, true);
            app.WindowLengthEditField.Position = [133 94 111 22];
            app.WindowLengthEditField.Value = 1000;

            % Create EpochDropDownLabel
            app.EpochDropDownLabel = uilabel(app.Panel);
            app.EpochDropDownLabel.HorizontalAlignment = 'right';
            app.EpochDropDownLabel.Position = [6 52 40 22];
            app.EpochDropDownLabel.Text = 'Epoch';

            % Create EpochDropDown
            app.EpochDropDown = uidropdown(app.Panel);
            app.EpochDropDown.Items = {};
            app.EpochDropDown.ValueChangedFcn = createCallbackFcn(app, @EpochDropDownValueChanged, true);
            app.EpochDropDown.Position = [61 52 183 22];
            app.EpochDropDown.Value = {};

            % Create BoutMoreRipplesLabel
            app.BoutMoreRipplesLabel = uilabel(app.Panel);
            app.BoutMoreRipplesLabel.Position = [398 94 131 22];
            app.BoutMoreRipplesLabel.Text = 'Bout with more ripples: ';

            % Create TotalPFCripplesLabel
            app.TotalPFCripplesLabel = uilabel(app.Panel);
            app.TotalPFCripplesLabel.Position = [399 69 100 22];
            app.TotalPFCripplesLabel.Text = 'Total PFC ripples:';

            % Create TotalHPCripplesLabel
            app.TotalHPCripplesLabel = uilabel(app.Panel);
            app.TotalHPCripplesLabel.Position = [399 44 102 22];
            app.TotalHPCripplesLabel.Text = 'Total HPC ripples:';

            % Create BoutDropDownLabel
            app.BoutDropDownLabel = uilabel(app.Panel);
            app.BoutDropDownLabel.HorizontalAlignment = 'right';
            app.BoutDropDownLabel.Position = [6 10 30 22];
            app.BoutDropDownLabel.Text = 'Bout';

            % Create BoutDropDown
            app.BoutDropDown = uidropdown(app.Panel);
            app.BoutDropDown.Items = {};
            app.BoutDropDown.ValueChangedFcn = createCallbackFcn(app, @BoutDropDownValueChanged, true);
            app.BoutDropDown.Position = [51 10 67 22];
            app.BoutDropDown.Value = {};

            % Create EventDropDownLabel
            app.EventDropDownLabel = uilabel(app.Panel);
            app.EventDropDownLabel.HorizontalAlignment = 'right';
            app.EventDropDownLabel.Position = [127 10 36 22];
            app.EventDropDownLabel.Text = 'Event';

            % Create EventDropDown
            app.EventDropDown = uidropdown(app.Panel);
            app.EventDropDown.Items = {};
            app.EventDropDown.ValueChangedFcn = createCallbackFcn(app, @EventDropDownValueChanged, true);
            app.EventDropDown.Position = [178 10 66 22];
            app.EventDropDown.Value = {};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = RipplesAnalyzerMatFile

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end