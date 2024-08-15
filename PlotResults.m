function PlotResults (time, y_out, y_dbm, cleaned_dataset, WhichResult,path,SaveOption,GlucoseCtrlOptn,TreatmentTime,SimulationTime)
exp_data_switch = 1;
RunNumber = '';
colors = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840]];
days_to_wks = 1/7;
[AGE_data,MCP_data,MCP_error] = AGEMCPData ();
M_struct = GlucoseDataStructuring();
Total_FC_dbdb = (y_out)./y_dbm(end,:);

    % Glucose input
    if sum(strcmp(WhichResult,'Glucose'))>0
        figure(60)
        Weeks2Days = 7;
        SimTimeDays = SimulationTime*Weeks2Days;
        timespan = linspace(1,SimTimeDays,SimTimeDays);
        for i = 1:length(timespan)
            GlucoseConcentration(i) = GlucoseInputFunc(time(i),GlucoseCtrlOptn,TreatmentTime,SimulationTime);
        end
        plot(time*days_to_wks,GlucoseConcentration/GlucoseConcentration(1),'linewidth',3,'Color','black')
        hold on
        % Plot the glucose experimental data
        for i = 1:6
            errorbar(M_struct{i}(1,:),M_struct{i}(2,:),M_struct{i}(3,:),'o','Color','black', 'HandleVisibility','off')
            scatter(M_struct{i}(1,:),M_struct{i}(2,:),'filled')
            hold on
        end
        xlabel("Time (weeks)"), ylabel("Glucose (g/ml)"),grid on
        legend('Model input','Cohen2001','Kolavennu2006','Ziyadeh2000','Koya2000','Cohen1995','Cohen2002'),legend('Location','best')
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        full_path_glucose = strcat(path,'\glucose_data_fit');
        export_fig(full_path_glucose,'-r1000','-a4', '-q101', '-painters', '-png')
    end
    
    % AGE
    if sum(strcmp(WhichResult,'AGE'))>0  
        species = 1;
        figure(70)
        plot(time*days_to_wks,Total_FC_dbdb(:,1),'linewidth',3,'LineStyle','-.','Color','black'),xlabel('Time (weeks)'),ylabel('AGE fold Change'),grid on, ylim([0 1.25*max(Total_FC_dbdb(:,1))])
        hold on
        scatter(cleaned_dataset {species}(:,1)*days_to_wks,cleaned_dataset {species}(:,4),'filled','MarkerFaceColor','black')
        errorbar(cleaned_dataset{species}(:,1)*days_to_wks,cleaned_dataset{species}(:,4),cleaned_dataset{species}(:,5),'o','MarkerEdgeColor','black', 'Color', 'black')
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial'); legend('Model','Exp data')
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_age_fc = strcat(path,'\AGE_fc');
        export_fig(full_path_age_fc,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
    
    % MCP
    if sum(strcmp(WhichResult,'MCP'))>0   
        figure (80)
        plot(time*days_to_wks,Total_FC_dbdb(:,2),'linewidth',3,'Color',colors(1,:)),xlabel('Time (weeks)'),ylabel('MCP-1 fold change'),grid on, ylim([0 1.25*max(Total_FC_dbdb(:,2))])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_mcp_fc = strcat(path,'\MCPfc');
        export_fig(full_path_mcp_fc,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end

% Macrophage fold change
if sum(strcmp(WhichResult,'Macrophage'))>0 || sum(strcmp(WhichResult,'All'))>0
    figure(10)
    if exp_data_switch == 1
        species_M = 2;
        Varybl = 3;
        plot(time*days_to_wks,Total_FC_dbdb(:,3),'linewidth',3,'Color',colors(2,:)),xlabel('Time (weeks)'),ylabel('Macrophage fold change'),grid on, ylim([0 1.25*max(Total_FC_dbdb(:,Varybl))])
        hold on
        scatter(cleaned_dataset {species_M}(:,1)*days_to_wks,cleaned_dataset {species_M}(:,4),'filled','MarkerFaceColor',colors(2,:))
        errorbar(cleaned_dataset{species_M}(:,1)*days_to_wks,cleaned_dataset{species_M}(:,4),cleaned_dataset{species_M}(:,5),'o','MarkerEdgeColor',colors(2,:), 'Color', colors(2,:))
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');legend('Model','Exp data')
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_macrophage_fc = strcat(path,'\Macrophage_fc',RunNumber);
        export_fig(full_path_macrophage_fc,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
end

% TGF-B fold change
if sum(strcmp(WhichResult,'TGF'))>0 || sum(strcmp(WhichResult,'All'))>0
    figure(20)
    if exp_data_switch == 1
        Varybl = 4;
        species_M = 3;
        plot(time*days_to_wks,Total_FC_dbdb(:,Varybl),'linewidth',3,'Color',colors(3,:)),xlabel('Time (weeks)'),ylabel('TGF fold Change'),grid on, ylim([0 1.25*max(Total_FC_dbdb(:,Varybl))])
        hold on
        scatter(cleaned_dataset {species_M}(:,1)*days_to_wks,cleaned_dataset {species_M}(:,4),'filled','MarkerFaceColor',colors(3,:))
        errorbar(cleaned_dataset{species_M}(:,1)*days_to_wks,cleaned_dataset{species_M}(:,4),cleaned_dataset{species_M}(:,5),'o','MarkerEdgeColor',colors(3,:), 'Color', colors(3,:))
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial'),legend('Model','Exp data')
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_tgf_fc = strcat(path,'\tgf_fc',RunNumber);
        export_fig(full_path_tgf_fc,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
end

% Activated mesangial cells fold change
if sum(strcmp(WhichResult,'AMC'))>0 || sum(strcmp(WhichResult,'All'))>0
    figure(30)
    if exp_data_switch == 1
        Varybl = 5;
        species_ma = 4;
        plot(time*days_to_wks,Total_FC_dbdb(:,Varybl),'linewidth',3,'Color',colors(4,:)),grid on
        ylim([0 1.25*max(Total_FC_dbdb(:,Varybl))])
        hold on
        scatter(cleaned_dataset{species_ma}(:,1)*days_to_wks,cleaned_dataset{species_ma}(:,4),'filled','MarkerFaceColor',colors(4,:))
        errorbar(cleaned_dataset{species_ma}(:,1)*days_to_wks,cleaned_dataset{species_ma}(:,4),cleaned_dataset{species_ma}(:,5),'o','MarkerEdgeColor',colors(4,:), 'Color', colors(4,:))
        xlabel("Time (weeks)"), ylabel("AMC fold change"),grid on,legend('Model','Exp data'),legend('Location','best'),%title("Activated Mesangial")
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_mafc = strcat(path,'\AMC_fc');
        export_fig(full_path_mafc,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
end
% Collagen fold change
if sum(strcmp(WhichResult,'Collagen'))>0 || sum(strcmp(WhichResult,'All'))>0
    figure(40)
    if exp_data_switch == 1
        Varybl = 8;
        species_rho = 5;
        plot(time*days_to_wks,Total_FC_dbdb(:,Varybl),'linewidth',3,'Color',colors(6,:)),grid on,ylim([0 1.25*max(Total_FC_dbdb(:,Varybl))])
        hold on
        scatter(cleaned_dataset{species_rho}(:,1)*days_to_wks,cleaned_dataset{species_rho}(:,4),'filled','MarkerFaceColor',colors(6,:))
        errorbar(cleaned_dataset{species_rho}(:,1)*days_to_wks,cleaned_dataset{species_rho}(:,4),cleaned_dataset{species_rho}(:,5),'o','MarkerEdgeColor',colors(6,:), 'Color', colors(6,:))
        xlabel("Time (weeks)"), ylabel("Collagen fold change"),legend('Model','Exp data'),grid on,legend('Location','best'),%title("Collagen")
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_rhofc = strcat(path,'\Collagen_fc');
        export_fig(full_path_rhofc,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
end   

% MMP
if sum(strcmp(WhichResult,'MMP'))>0 
    figure (80)
    Varybl = 6;
    plot(time*days_to_wks,Total_FC_dbdb(:,Varybl),'linewidth',3,'Color',colors(5,:)),grid on
    xlabel("Time (weeks)"), ylabel("MMP fold change"),ylim([0 1.25*max(Total_FC_dbdb(:,Varybl))])
    get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
    if sum(strcmp(SaveOption,'Save'))>0
    full_path_mmpfc = strcat(path,'\MMPFC');
    export_fig(full_path_mmpfc,'-r1000','-a4', '-q101', '-painters', '-png')
    end
end
% TIMP
if sum(strcmp(WhichResult,'TIMP'))>0 
    figure (90)
    Varybl = 7;
    plot(time*days_to_wks,Total_FC_dbdb(:,Varybl),'linewidth',3,'Color',colors(7,:)),grid on
    xlabel("Time (weeks)"), ylabel("TIMP fold change"),ylim([0 1.25*max(Total_FC_dbdb(:,Varybl))])
    get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
    if sum(strcmp(SaveOption,'Save'))>0
    full_path_timpfc = strcat(path,'\TIMPfc');
    export_fig(full_path_timpfc,'-r1000','-a4', '-q101', '-painters', '-png')
    end
end

% AGE-MCP dynamics with experimental data
if sum(strcmp(WhichResult,'AGEMCP'))>0
    if exp_data_switch == 1
        figure(50)
        plot(y_out(:,1),y_out(:,2),'linewidth',3,'Color',colors(1,:))
        hold on
        scatter(AGE_data,MCP_data ,'filled','MarkerFaceColor',colors(1,:)),grid on %ylim([0 3e-5])
        errorbar(AGE_data,MCP_data,MCP_error,'o','MarkerEdgeColor',colors(1,:), 'Color', colors(1,:))
        xlabel('AGE (g/ml)'),ylabel('MCP-1 (g/ml)'),legend('Model','Exp data'),get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_agemcp = strcat(path,'\AGE_MCP');
        export_fig(full_path_agemcp,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
end

% Absolute values plotted
    % Glucose input
    if sum(strcmp(WhichResult,'AbsGlucose'))>0
        Weeks2Days = 7;
        SimTimeDays = SimulationTime*Weeks2Days;
        timespan = linspace(1,SimTimeDays,SimTimeDays);
        for i = 1:length(timespan)
            GlucoseConcentration(i) = GlucoseInputFunc(timespan(i),GlucoseCtrlOptn,TreatmentTime);
        end
        
        figure(1)
        plot(timespan*days_to_wks,GlucoseConcentration,'linewidth',3,'Color','black')
        ylabel("Glucose (g/ml)"),xlabel("Time (weeks)"),grid on, ylim([0 1.25*max(GlucoseConcentration)])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_glucose = strcat(path,'\Glucose');
        export_fig(full_path_glucose,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
    % AGE
    if sum(strcmp(WhichResult,'AbsAGE'))>0   
        figure(2)
        plot(time*days_to_wks,y_out(:,1),'linewidth',3,'LineStyle','-.','Color','black'),xlabel('Time (weeks)'),ylabel('AGE (g/ml)'),grid on, ylim([0 1.25*max(y_out(:,1))]) %axis([10 25 0 1]),
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_age = strcat(path,'\AGE');
        export_fig(full_path_age,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
    % MCP
    if sum(strcmp(WhichResult,'AbsMCP'))>0   
        figure (3)
        plot(time*days_to_wks,y_out(:,2),'linewidth',3,'Color',colors(1,:)),xlabel('Time (weeks)'),ylabel('MCP-1 (g/ml)'),grid on, ylim([0 1.25*max(y_out(:,2))]) %axis([10 25 0.002 0.0025]),
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_mcp = strcat(path,'\MCP');
        export_fig(full_path_mcp,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
    % Macrophage
    if sum(strcmp(WhichResult,'AbsMacrophage'))>0
        figure(4)
        plot(time*days_to_wks,y_out(:,3),'linewidth',3,'Color',colors(2,:)),grid on
        xlabel("Time (weeks)"), ylabel("Macrophage (g/ml)"),ylim([0 1.25*max(y_out(:,3))])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_macrophage = strcat(path,'\Macrophage');
        export_fig(full_path_macrophage,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
    % TGF-B
    if sum(strcmp(WhichResult,'AbsTGF'))>0    
        figure (5)
        plot(time*days_to_wks,y_out(:,4),'linewidth',3,'Color',colors(3,:)),grid on
        xlabel("Time (weeks)"), ylabel("TGF-B (g/ml)"), ylim([0 1.25*max(y_out(:,4))])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_tgfb = strcat(path,'\TGF-B');
        export_fig(full_path_tgfb,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end

    % AMC
    if sum(strcmp(WhichResult,'AbsAMC'))>0 
        figure(6)
        plot(time*days_to_wks,y_out(:,5),'linewidth',3,'Color',colors(4,:)),grid on
        xlabel("Time (weeks)"), ylabel("AMC (g/ml)"),ylim([0 1.25*max(y_out(:,5))])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_ma = strcat(path,'\AMC');
        export_fig(full_path_ma,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end

    % MMP
    if sum(strcmp(WhichResult,'AbsMMP'))>0 
        figure (7)
        plot(time*days_to_wks,y_out(:,6),'linewidth',3,'Color',colors(5,:)),grid on
        xlabel("Time (weeks)"), ylabel("MMP (g/ml)"),ylim([0 1.25*max(y_out(:,6))])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_mmp = strcat(path,'\MMP');
        export_fig(full_path_mmp,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
    % TIMP
    if sum(strcmp(WhichResult,'AbsTIMP'))>0 
        figure (8)
        plot(time*days_to_wks,y_out(:,7),'linewidth',3,'Color',colors(7,:)),grid on
        xlabel("Time (weeks)"), ylabel("TIMP (g/ml)"),ylim([0 1.25*max(y_out(:,7))])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_timp = strcat(path,'\TIMP');
        export_fig(full_path_timp,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end

    % collagen
    if sum(strcmp(WhichResult,'AbsCollagen'))>0 
        figure (9)
        plot(time*days_to_wks,y_out(:,8),'linewidth',3,'Color',colors(6,:)),grid on
        xlabel("Time (weeks)"), ylabel("Collagen (g/ml)"),ylim([0 1.25*max(y_out(:,8))])
        get(gca);set(gcf, 'Color', 'w');set(gca,'FontSize',14,'FontName','Arial');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_rho = strcat(path,'\Collagen');
        export_fig(full_path_rho,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end

    % Plot remainder species
    if sum(strcmp(WhichResult,'RemainderSpeciesSubPlot'))>0   
        linestyle = {'-.','-','-'};
        y_plot = Total_FC_dbdb;
        timespan_dbm = linspace(1,168,168);
        figure(100)
        subplot(3,1,1)
        plot(timespan_dbm*days_to_wks,GlucoseConcentration/GlucoseConcentration(1),'linewidth',1.5,'Color','k')
        ylabel('Fold change'),legend('Glucose'),%ylim([0 5])
        subplot(3,1,2)
        plot(time*days_to_wks,y_plot(:,1),'linewidth',1.5,'Color','k','LineStyle','-.')
        hold on
        for i = 2:3
            plot(time*days_to_wks,y_plot(:,i),'linewidth',1.5,'Color',colors(i-1,:))
            hold on
        end
        ylabel('Fold change'),xlabel('Time (weeks)'),%ylim([0 5])
        legend("AGE","MCP","MAC")
        subplot(3,1,3)
        for i = 6:7
            plot(time*days_to_wks,y_plot(:,i),'linewidth',1.5,'Color',colors(i-2,:))
            hold on 
        end
        ylabel('Fold change'),xlabel('Time (weeks)'),%ylim([0 5])
        legend("MMP","TIMP")
        get(gca);set(gcf, 'Color', 'w');
        if sum(strcmp(SaveOption,'Save'))>0
        full_path_var_fc = strcat(path,'\var_fc ',RunNumber);
        export_fig(full_path_var_fc ,'-r1000','-a4', '-q101', '-painters', '-png')
        end
    end
end
