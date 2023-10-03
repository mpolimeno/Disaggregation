% to compute relevant metrics

% we can make it smarter and avoid all of these loops now that we have the
% data

clear
close all
clc
SEED_vec = 1:50;
NC_vec = [10,25,50,100,150,200,250,300,400];
FLOW = 2;

Data_stresses = [];

NC_0_list  = [];
Rg0_list   = [];
Rmax0_list = [];


NC_1_list  = [];
Rg1_list   = [];
Rmax1_list = [];

NC_2_list  = [];
Rg2_list   = [];
Rmax2_list = [];

loop_counter = 1;

stress_distance = [];
stress_location = [];

Data_sizes = [];

for ii=1:length(NC_vec)
    NC = NC_vec(ii);
    stress_distr_NC = [];
    for jj=1:length(SEED_vec)
        SEED = SEED_vec(jj);
         if isfile(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_NC_%i_BROKEN_FLOW_%i.mat",SEED,NC,FLOW))==1
             
            % Here we get the position of the cubes in the original
            % aggregate
            F0 = load(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_NC_%i_BROKEN_FLOW_%i.mat",SEED,NC,FLOW));
            pos_0 = F0.xc;
            NC_0_list = [NC_0_list size(pos_0,1)];
            if size(pos_0,1)==1
                cm = pos_0; % for one cube only
            else
                cm = mean(pos_0); % for multiple cubes
            end
            % Store all the Radii of Gyration
            Rg_0 = compute_Rg(pos_0,cm);
            Rg0_list = [Rg0_list; Rg_0];
            
            % Store all of the Maximum Radii
            Rmax_0 = compute_Rmax(pos_0,cm);
            Rmax0_list = [Rmax0_list; Rmax_0];
            
            % Here we get which two cubes will break
            Cubes = load(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_NC_%i_BROKEN_FLOW_%i.mat",SEED,NC,FLOW));
            Cubes_to_break = Cubes.cubes_that_break;

            com_1 = pos_0(Cubes_to_break(1),:);
            com_2 = pos_0(Cubes_to_break(2),:);

            face_location = 0.5*(com_2-com_1);


            stress_dist = norm(cm-face_location);
            stress_distance = [stress_distance; NC Rg_0 Rmax_0 stress_dist];

            % Here we get the position of the cubes in the first aggregate
            F1 = load(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_NC_%i_BROKEN_FLOW_%i.mat",SEED,NC,FLOW));
            pos_1 = F1.xc_1;
            NC_1_list = [NC_1_list size(pos_1,1)];
            if size(pos_1,1)==1
                cm = pos_1; % for one cube only
            else
                cm = mean(pos_1); % for multiple cubes
            end
            % Store all the Radii of Gyration
            Rg_1 = compute_Rg(pos_1,cm);
            Rg1_list = [Rg1_list; Rg_1];
            
            % Store all of the Maximum Radii
            Rmax_1 = compute_Rmax(pos_1,cm);
            Rmax1_list = [Rmax1_list; Rmax_1];
            NC_1 = size(pos_1,1);

            % Here we get the position of the cubes in the second aggregate
            F2 = load(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_NC_%i_BROKEN_FLOW_%i.mat",SEED,NC,FLOW));
            pos_2 = F2.xc_2;
            NC_2_list = [NC_2_list size(pos_2,1)];
            if size(pos_2,1)==1
                cm = pos_2; % for one cube only
            else
                cm = mean(pos_2); % for multiple cubes
            end
            % Store all the Radii of Gyration
            Rg_2 = compute_Rg(pos_2,cm);
            Rg2_list = [Rg2_list; Rg_2];
            
            % Store all of the Maximum Radii
            Rmax_2 = compute_Rmax(pos_2,cm);
            Rmax2_list = [Rmax2_list; Rmax_2];
            
            NC_2 = size(pos_2,1);
        end
        
        if isfile(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_Internal_stresses_NC_%i_FLOW_%i.mat",SEED,NC,FLOW))==1            
            F_stress = load(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_Internal_stresses_NC_%i_FLOW_%i.mat",SEED,NC,FLOW));
            
            internal_stresses = F_stress.internal_stresses;
            norms_list = [];
            for mm=1:size(internal_stresses,2)
                norm_of_stresses = norm(internal_stresses(:,mm));
                norms_list = [norms_list;norm_of_stresses];
            end
            max_internal_stress = max(norms_list);
            mean_internal_stress = mean(norms_list);
            Data_stresses = [Data_stresses; NC Rg_0 max_internal_stress mean_internal_stress Rmax_0];
            
            % I think it would be a good idea to store the aggregate-size
            % data in the following way:
            % Original aggregate data, larger of the two data, smaller of
            % the two data
            
            if (NC_1>=NC_2)
                Data_sizes = [Data_sizes; NC Rg_0 Rmax_0 NC_1 Rg_1 Rmax_1 NC_2 Rg_2 Rmax_2]; % This will be used for the comparison of sizes
            else
                Data_sizes = [Data_sizes; NC Rg_0 Rmax_0 NC_2 Rg_2 Rmax_2 NC_1 Rg_1 Rmax_1]; % This will be used for the comparison of sizes
            end
            
        elseif isfile(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_Internal_stresses_NC_%i_FLOW_%i_LoopEvent_%i.mat",SEED,NC,FLOW,loop_counter))==1
            F_stress = load(sprintf("~/Desktop/NEGOTIUM/UCM/Disaggregation/Codes_FB/IAA_2023_10/Rotation/Stats/SEED_%i_ORIGINAL_Agg_Internal_stresses_NC_%i_FLOW_%i_LoopEvent_%i.mat",SEED,NC,FLOW,loop_counter));

            % Will need to think about what we want to do about loop
            % events. But if they only happen rarely, I do not think
            % processing their position data is at all worth it as it would
            % mess up the sizes of the arrays
            pos_0 = F_stress.xc;
            NC_0_list = [NC_0_list size(pos_0,1)];
            if size(pos_0,1)==1
                cm = pos_0; % for one cube only
            else
                cm = mean(pos_0); % for multiple cubes
            end
            % Store all the Radii of Gyration
            Rg_0 = compute_Rg(pos_0,cm);
            Rg0_list = [Rg0_list; Rg_0];
            
            % Store all of the Maximum Radii
            Rmax_0 = compute_Rmax(pos_0,cm);
            Rmax0_list = [Rmax0_list; Rmax_0];
            
            loop_counter = loop_counter + 1;
            
            internal_stresses = F_stress.internal_stresses;
            norms_list = [];
            for mm=1:size(internal_stresses,2)
                norm_of_stresses = norm(internal_stresses(:,mm));
                norms_list = [norms_list;norm_of_stresses];
            end
            max_internal_stress = max(norms_list);
            mean_internal_stress = mean(norms_list);
            Data_stresses = [Data_stresses; NC Rg_0 max_internal_stress mean_internal_stress Rmax_0];
            
            % To have the arrays being the same sizes, we set the appended
            % values to 0, as the original aggregate is not broken
            NC_1_list = [NC_1_list 0];
            Rg1_list = [Rg1_list; 0];
            Rmax1_list = [Rmax1_list; 0];
            
            NC_2_list = [NC_2_list 0];
            Rg2_list = [Rg2_list; 0];
            Rmax2_list = [Rmax2_list; 0];
            
            Data_sizes = [Data_sizes; NC Rg_0 Rmax_0 0 0 0 0 0 0]; % To have consistent array sizes
        end
%         x = ones(size(norms_list,1),1)*Rg_0;
%         figure(25)
%         plot(x,norms_list,'.')
%         hold on
%         xlabel("$R_g$",'interpreter','latex')
%         ylabel("||f||")
    end
    stress_distr_NC = [stress_distr_NC,norms_list];
    figure(27+ii)
    histogram(stress_distr_NC)
end


% Because of loop events, I either need to add a zero to the NC_1 and NC_2
% arrays, but then make sure that M/M_1 or M/M_2 is set to 0 to avoid
% numerical issues. But the easier solution is to simply dismiss that data
% to begin with, to be honest

% Array with both Rg and correponding M
Data_Rg0_list = zeros(size(Rg0_list,1),2);
Data_Rg0_list(:,1) = Rg0_list;
Data_Rg0_list(:,2) = NC_0_list;

Rg_avg = [];
counter = 0;
for kk=1:size(Data_Rg0_list,1)
    if kk>1 && mod(kk,50)==0
        old_counter = counter+1;
        Rg_avg = [Rg_avg;mean(Data_Rg0_list(old_counter:kk,1))];
        counter = kk;
    end
end


% Array with both Rmax and corresponding M
Data_Rmax0_list = zeros(size(Rmax0_list,1),2);
Data_Rmax0_list(:,1) = Rmax0_list;
Data_Rmax0_list(:,2) = NC_0_list;

% Array with both Rg and correponding M
Data_Rg1_list = zeros(size(Rg1_list,1),2);
Data_Rg1_list(:,1) = Rg1_list;
Data_Rg1_list(:,2) = NC_1_list;

% Array with both Rmax and corresponding M
Data_Rmax1_list = zeros(size(Rmax1_list,1),2);
Data_Rmax1_list(:,1) = Rmax1_list;
Data_Rmax1_list(:,2) = NC_1_list;

% Array with both Rg and correponding M
Data_Rg2_list = zeros(size(Rg2_list,1),2);
Data_Rg2_list(:,1) = Rg2_list;
Data_Rg2_list(:,2) = NC_2_list;

% Array with both Rmax and corresponding M
Data_Rmax2_list = zeros(size(Rmax2_list,1),2);
Data_Rmax2_list(:,1) = Rmax2_list;
Data_Rmax2_list(:,2) = NC_2_list;

% For ratios (will clean up later, as things are redundant now)
% % if (NC_1>=NC_2)
% %     Data_sizes = [Data_sizes; NC Rg_0 Rmax_0 NC_1 Rg_1 Rmax_1 NC_2 Rg_2 Rmax_2]; % This will be used for the comparison of sizes
% % else
% %     Data_sizes = [Data_sizes; NC Rg_0 Rmax_0 NC_2 Rg_2 Rmax_2 NC_1 Rg_1 Rmax_1]; % This will be used for the comparison of sizes
% % end
M1_over_M   = zeros(size(Data_sizes,1),1);
M2_over_M   = zeros(size(Data_sizes,1),1);
Rg1_over_Rg = zeros(size(Data_sizes,1),1);
Rg2_over_Rg = zeros(size(Data_sizes,1),1);
for ii=1:size(Data_sizes,1)
    % Ratio between M of original aggregate, and larger of the two broken up
    M1_over_M(ii) = Data_sizes(ii,4)/Data_sizes(ii,1);

    % Ratio between M of original aggregate, and smaller of the two broken up
    M2_over_M(ii) = Data_sizes(ii,7)/Data_sizes(ii,1);
    
    % Ratio between Rg of original aggregate, and larger of the two broken up
    Rg1_over_Rg(ii) = Data_sizes(ii,5)/Data_sizes(ii,2);
   
    % Ratio between Rg of original aggregate, and smaller of the two broken up
    Rg2_over_Rg(ii) = Data_sizes(ii,8)/Data_sizes(ii,2);
end



%%% This one is too noisy; It would be better to have one plot per each M
figure(1)
plot(M1_over_M,Rg1_over_Rg,'r.','markersize',10)
hold on
plot(M2_over_M,Rg2_over_Rg,'k.','markersize',10)
legend("$R_{g1}/R_g$","$R_{g2}/R_g$",'interpreter','latex','location','southeast','fontsize',20)
ax = gca;
ax.FontSize = 25;
xlabel("$M_{s}/M$",'interpreter','latex')
ylabel("$R_{gs}/R_g$",'interpreter','latex')
title("Total Relative Sizes",'interpreter','latex','fontsize',20)

% % I think adding also the original aggregates would be useful for
% % comparison purposes
figure(2)
loglog(Data_Rg0_list(:,2),Data_Rg0_list(:,1),'m.','markersize',10)
ax = gca;
ax.FontSize = 25;
xlabel("$M$",'interpreter','latex')
ylabel("$R_g$",'interpreter','latex')
xlim([1 400])
ylim([0 20])
title("Radius of Gyration",'interpreter','latex','fontsize',20)

figure(3)
loglog(Data_Rg1_list(:,2),Data_Rg1_list(:,1),'m.','markersize',10)
ax = gca;
ax.FontSize = 25;
xlabel("$M_1$",'interpreter','latex')
ylabel("$R_{g1}$",'interpreter','latex')
xlim([1 400])
ylim([0 20])
title("Radius of Gyration",'interpreter','latex','fontsize',20)


figure(4)
loglog(Data_Rg2_list(:,2),Data_Rg2_list(:,1),'m.','markersize',10)
ax = gca;
ax.FontSize = 25;
xlabel("$M_2$",'interpreter','latex')
ylabel("$R_{g2}$",'interpreter','latex')
xlim([1 400])
ylim([0 20])
title("Radius of Gyration",'interpreter','latex','fontsize',20)

%%%% MAX RADIUS (for later, if needed)
% 
% figure(3)
% loglog(Data_Rmax0_list(:,2),Data_Rmax0_list(:,1),'m.','markersize',10)
% hold on
% loglog(Data_Rmax1_list(:,2),Data_Rmax1_list(:,1),'k.','markersize',10)
% hold on
% loglog(Data_Rmax2_list(:,2),Data_Rmax2_list(:,1),'r.','markersize',10)
% legend("Original Aggregate","Aggregate 1","Aggregate 2",'location','southeast','fontsize',20)
% ax = gca;
% ax.FontSize = 25;
% xlabel("$M$",'interpreter','latex')
% ylabel("$R_{max}$",'interpreter','latex')
% title("Maximum Radius",'interpreter','latex','fontsize',20)

%%% MAX STRESS NO FIT
% 
% figure(3)
% plot(Data_stresses(:,2),Data_stresses(:,3),'b.','markersize',10)
% ax = gca;
% ax.FontSize = 25;
% xlabel("$R_g$",'interpreter','latex')
% ylabel("$\max{||\vec{f}_{\mathrm {internal}}||}$",'interpreter','latex')
% title("Maximum Internal Stress",'interpreter','latex','fontsize',20)
%


figure(5)
plot(Data_stresses(:,2),Data_stresses(:,4),'m.','markersize',10)
ax = gca;
ax.FontSize = 25;
xlabel("$R_g$",'interpreter','latex')
ylabel("$mean{||\vec{f}_{\mathrm {internal}}||}$",'interpreter','latex')
title("Mean Internal Stress (Linear Scale)",'interpreter','latex','fontsize',20)

figure(23)
plot(Data_stresses(:,1),Data_stresses(:,4),'m.','markersize',10)
ax = gca;
ax.FontSize = 25;
xlabel("$M$",'interpreter','latex')
ylabel("$mean{||\vec{f}_{\mathrm {internal}}||}$",'interpreter','latex')
title("Mean Internal Stress (Linear Scale)",'interpreter','latex','fontsize',20)

cubic_fit = polyfit(Data_stresses(:,2),Data_stresses(:,3),3);

best_fit = cubic_fit(1)*Data_stresses(:,2).^3 + cubic_fit(2)*Data_stresses(:,2).^2 + cubic_fit(3)*Data_stresses(:,2) + cubic_fit(1);

y_cubed = cubic_fit(1)*Data_stresses(:,2).^3;

figure(6)
plot(Data_stresses(:,2),Data_stresses(:,3),'b.','markersize',10)
hold on
plot(Data_stresses(:,2),best_fit,'r-','linewidth',4)
xlim([min(Data_stresses(:,2)),max(Data_stresses(:,2))])
ax = gca;
ax.FontSize = 25;
%xlabel("$R_g$",'interpreter','latex')
%ylabel("$\max{||\vec{f}_{\mathrm {internal}}||}$",'interpreter','latex')
%title("Maximum Internal Stress (Linear Scale)",'interpreter','latex','fontsize',20)

log_fit = polyfit(log(Data_stresses(:,2)),log(Data_stresses(:,3)),1);

best_log_fit = exp(log_fit(2)).*Data_stresses(:,2).^log_fit(1);

log_fit_M = polyfit(log(Data_stresses(:,1)),log(Data_stresses(:,3)),1);

best_log_fit_M = exp(log_fit_M(2)).*Data_stresses(:,1).^log_fit_M(1);

figure(7)
loglog(Data_stresses(:,2),Data_stresses(:,3),'b.','markersize',10)
hold on
loglog(Data_stresses(:,2),best_log_fit,'r-','linewidth',4)
hold on
loglog(Data_stresses(:,2),y_cubed,'k--','linewidth',2)
ax = gca;
ax.FontSize = 25;
%xlabel("$R_g$",'interpreter','latex')
%ylabel("$\max{||\vec{f}_{\mathrm {internal}}||}$",'interpreter','latex')
%title("Maximum Internal Stress (Log-Log scale)",'interpreter','latex','fontsize',20)


figure(22)
loglog(Data_stresses(:,1),Data_stresses(:,3),'b.','markersize',10)
hold on
loglog(Data_stresses(:,1),best_log_fit_M,'r-','linewidth',4)
ax = gca;
ax.FontSize = 25;
xlabel("$M$",'interpreter','latex')
ylabel("$\max{||\vec{f}_{\mathrm {internal}}||}$",'interpreter','latex')
title("Maximum Internal Stress (Log-Log scale)",'interpreter','latex','fontsize',20)

counter = 1;
index_list = zeros(length(NC_vec)-1,1);
while(counter<length(NC_vec))
    for kk=1:size(stress_distance,1)
        if stress_distance(kk+1)==stress_distance(kk)
            M_value = stress_distance(kk,1);
        else
            index_to_end = kk;
            index_list(counter) = index_to_end;
            figure(7+counter)
            if counter==1
                histogram(stress_distance(1:index_to_end,4),"Normalization","PDF")
                ax = gca;
                ax.FontSize = 25;
                xlabel("d","interpreter","latex")
                ylabel("Relative Frequency","interpreter","latex")
                title(sprintf("Distribution of Distances for M=%i",M_value),'interpreter','latex')
            else
                histogram(stress_distance(index_list(counter-1)+1:index_to_end,4),"Normalization","PDF")
                ax = gca;
                ax.FontSize = 25;
                xlabel("d","interpreter","latex")
                ylabel("Relative Frequency","interpreter","latex")
                title(sprintf("Distribution of Distances for M=%i",M_value),'interpreter','latex')
            end
            counter = counter + 1;
        end
    end
end

counter_2 = 1;
ending = length(SEED_vec)*length(NC_vec);
for nn=1:ending
    figure(counter+7+counter_2)
    if mod(nn,length(SEED_vec))==0
        M_original = NC_vec(counter_2);
        if nn==length(SEED_vec)
            plot(M1_over_M(1:nn),Rg1_over_Rg(1:nn),'r.','markersize',10)
            hold on
            plot(M2_over_M(1:nn),Rg2_over_Rg(1:nn),'k.','markersize',10)
            legend("$R_{g1}/R_g$","$R_{g2}/R_g$",'interpreter','latex','location','southeast','fontsize',20)
            ax = gca;
            ax.FontSize = 25;
            xlabel("$M_{s}/M$",'interpreter','latex')
            ylabel("$R_{gs}/R_g$",'interpreter','latex')
            title(sprintf("Relative Sizes for M=%i",M_original),'interpreter','latex')
        else
            plot(M1_over_M(old_nn+1:nn),Rg1_over_Rg(old_nn+1:nn),'r.','markersize',10)
            hold on
            plot(M2_over_M(old_nn+1:nn),Rg2_over_Rg(old_nn+1:nn),'k.','markersize',10)
            legend("$R_{g1}/R_g$","$R_{g2}/R_g$",'interpreter','latex','location','southeast','fontsize',20)
            ax = gca;
            ax.FontSize = 25;
            xlabel("$M_{s}/M$",'interpreter','latex')
            ylabel("$R_{gs}/R_g$",'interpreter','latex')
            title(sprintf("Relative Sizes for M=%i",M_original),'interpreter','latex')
        end
        old_nn = nn;
        counter_2 = counter_2 + 1;
    end
end


counter_3 = 1;
ending = length(SEED_vec)*length(NC_vec);
%cmap = hot;
%color_list = ["#0072BD","#FF0000","#FF00FF","#0000FF","#00FFFF","#77AC30","#000000"];
color_list = ["#0072BD","#FF0000","#FF00FF","#0000FF","#00FFFF","#77AC30","#000000","#D95319","#EDB120"];
color_count = 0;
leg = [];
for mm=1:ending
    figure(counter+7+counter_2)
    if mod(mm,length(SEED_vec))==0
        %index = counter_3*20;
        %col = cmap(index,:);
        color_count = color_count+1;
        col = color_list(color_count);
        if mm==length(SEED_vec)
            plot(M1_over_M(1:mm),Rg1_over_Rg(1:mm),'.','color',col,'markersize',10)
            hold on
            plot(M2_over_M(1:mm),Rg2_over_Rg(1:mm),'.','color',col,'markersize',10)
            legend("$R_{g1}/R_g$","$R_{g2}/R_g$",'interpreter','latex','location','southeast','fontsize',20)
            legend off
            ax = gca;
            ax.FontSize = 25;
            xlabel("$M_{s}/M$",'interpreter','latex')
            ylabel("$R_{gs}/R_g$",'interpreter','latex')
            %legend_inputs = sprintf("M=%i",NC_vec(color_count));
            %leg = [leg legend_inputs];
            %legend(leg)
        else
            plot(M1_over_M(old_mm+1:mm),Rg1_over_Rg(old_mm+1:mm),'.','color',col,'markersize',10)
            hold on
            plot(M2_over_M(old_mm+1:mm),Rg2_over_Rg(old_mm+1:mm),'.','color',col,'markersize',10)
            legend("$R_{g1}/R_g$","$R_{g2}/R_g$",'interpreter','latex','location','southeast','fontsize',20)
            legend off
            ax = gca;
            ax.FontSize = 25;
            xlabel("$M_{s}/M$",'interpreter','latex')
            ylabel("$R_{gs}/R_g$",'interpreter','latex')
            %legend_inputs = sprintf("M=%i",NC_vec(color_count));
            %leg = [leg legend_inputs];
            %legend(leg)
        end
        old_mm = mm;
        counter_3 = counter_3 + 1;
    end
end
hold on
fractal_dim_IAA = 2.3;
M_ratio_vec = linspace(0,1,size(M1_over_M,1));
y = M_ratio_vec.^(1/2.2);
plot(M_ratio_vec,y,'linewidth',3)
% linear_fit_Rmax = polyfit(Data_stresses(:,5),Data_stresses(:,3),1);
% 
% best_fit_Rmax = linear_fit_Rmax(1)*Data_stresses(:,5) + linear_fit_Rmax(2);
% % 
% figure(8)
% plot(Data_stresses(:,5),Data_stresses(:,3),'r.','markersize',10)
% hold on
% plot(Data_stresses(:,5),best_fit_Rmax,'k-','linewidth',4)
% ax = gca;
% ax.FontSize = 25;
% xlabel("$R_{max}$",'interpreter','latex')
% ylabel("$\max{||\vec{f}_{\mathrm {internal}}||}$",'interpreter','latex')
% title("Maximum Internal Stress",'interpreter','latex','fontsize',20)

% Rg   = compute_Rg(aggregate,cm);
% Rmax = compute_Rmax(aggregate,cm);

% Radius of Gyration
function Rg = compute_Rg(aggregate,cm)
    [NC,dim] = size(aggregate);
    com = cm;
    tmp = 0;
    for ii = 1:NC
        for d=1:dim
            tmp = tmp + ((aggregate(ii,d)-com(d)))^2;
        end 
    end
    Rg = sqrt(tmp/NC);
end

function Rmax = compute_Rmax(aggregate,cm)
    [NC,dim] = size(aggregate);
    com = cm;
    tmp = zeros(NC,3);
    dist_max = zeros(NC,1);
    for ii = 1:NC
        for d=1:dim
            tmp(ii,d) = aggregate(ii,d)-com(d);
        end
        dist_max(ii) = norm(tmp(ii,:));
    end
    Rmax = 1 + max(dist_max);
end