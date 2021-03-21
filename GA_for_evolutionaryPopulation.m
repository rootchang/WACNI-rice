%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GA_for_evolutionaryPopulation.m: This script constructs a genetic algorithm to generate an evolutionary population targeting maximum of grain yield by 
%                                  changing values of 32 model parameters (28 biochemical/biophysical parameters and 4 plant size related parameters)
%
%  Input:
%   -
%
%  Output:
%   genotype_evolution_GW_GN.txt: file recording parameter values, grain weight and grain nitrogen concentration of all individuals
%   genotype_best.txt: file recording parameter values, grain weight and grain nitrogen concentration of the best individual in each generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxPoolNum=8; % CPU cores used for multi-processing
% maxNumCompThreads(maxPoolNum);
rng(sum(clock)*100);

% paprameters for genetic algorithm
population_size=100; %100
generation=100; %50
parents_num=10; %10
mutation_rate=5; %5
 
% size of plant organs
Root_DW=0.5; % g/tiller
Grain_NUM=160; % /ear
Grain_VOL=20.267; % mm3/grain
Leaf_DW=1; % g/tiller
stem_starch=0.75; % g/tiller
stem_maxVolume=1.8; % g/tiller
Area_leaf=0.0200; % m2/tiller
Thick_leaf=0.000175; % m
Vol_Fresh_ratio_leaf=0.25;
conc_organic_N_leaf=100; % mol/m3
conc_inorganic_N_leaf=50; % mol/m3
SLA0=0.02; % m2/g
leaf_structural_N=1/200; % g/g
conc_leaf_pro=0.035*6.25; % g/g
conc_inorganic_N_root=20; % mol/m3
conc_organic_N_root=10; % mol/m3
root_structural_N=1/200; % g/g
Vol_Fresh_ratio_root=0.25;
Fresh_dry_ratio_root=0.1;
Stem_pro=0.02; % g

PARA_list=[...
           [1,1];[1,1];[0.5,2];[1,1];[0.5,2];... % invariant: Vc_J, PPFD, Phi_CO2
           [1,1];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: lower_car_photo
           [0.5,2];[0.5,2];[0.5,2];[1,1];[0.5,2];... % invariant: conc_lower_organic_N_to_leaf_protein
           [1,1];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: gamma_star, 
           [0.5,2];[0.5,2];[1,1];[1,1];[1,1];... % invariant: grain_width_max, grain_length_width_ratio, Starch_stem
           [0.5,2];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: 
           [1,1];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: conc_inorganic_N_soil
           [0.5,2];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: 
           [0.5,2];[1,1];[1,1];[1,1];[1,1];[1,1];[1,1]...% invariant: tillerNumber; CO2a; protein_content_leaf; grain_structural_N; stem_structural_mass; Protein_stem
          ]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Names of parameters in PARA_list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Line1(01-05)  coef_Vc_J; coef_PPFD; coef_Area_leaf; coef_Phi_CO2; coef_K_protein_organic_N_leaf_cat (protein degradation); ...
%Line2(06-10)  coef_lower_car_photo; coef_K_organic_N_leaf_loading_cat; coef_K_suc_leaf_loading_cat; coef_K_N_assimilation_leaf_cat; coef_K_TP_suc_leaf_cat; ...
%Line3(11-15)  coef_K_TP_starch_leaf_cat; coef_K_starch_suc_leaf_cat; coef_K_organic_N_protein_leaf_cat; coef_conc_lower_organic_N_to_leaf_protein; coef_K_inorganic_N_leaf_unloading_cat; ...
%Line4(16-20)  coef_gamma_star; coef_grain_number; coef_K_grain_growth_cat; coef_K_grain_starch_store_cat; coef_K_suc_grain_unloading_cat; 
%Line5(21-25)  coef_K_grain_protein_store_cat; coef_K_organic_N_grain_unloading_cat; coef_grain_width_max; coef_grain_length_width_ratio; coef_Starch_stem; 
%Line6(26-30)  coef_K_suc_starch_stem_cat; coef_K_starch_suc_stem_cat; coef_K_organic_N_protein_stem_cat; coef_K_protein_organic_N_stem_cat; coef_Fresh_root; 
%Line7(31-35)  coef_conc_inorganic_N_soil; coef_K_inorganic_N_absorb_cat; coef_K_inorganic_N_root_loading_cat; coef_K_organic_N_root_loading_cat; coef_K_suc_root_unloading_cat; 
%Line8(36-40)  coef_K_organic_N_root_unloading_cat; coef_K_N_root_assimilation_cat; coef_K_root_growth_cat; coef_R_l_s_phloem0; coef_R_r_l_phloem0; 
%Line9(41-47)  coef_K_organic_N_xylem_to_lphloem_cat; coef_tillerNumber; coef_CO2a; protein_content_leaf; grain_structural_N; stem_structural_mass; Protein_stem

childen_genotype_list=[];
GW_sorted=[];
GN_sorted=[];
exist_childen_num=0;
if exist('genotype_evolution_GW_GN.txt', 'file')==2 % populations existing
    start_from_genotype=load('genotype_evolution_GW_GN.txt'); % GW, weight of all grains (grain yield); GN, grain nitrogen concentration
    exist_childen_num=size(start_from_genotype,1);
    disp(['Existing population size=',num2str(exist_childen_num)])
    if exist_childen_num>0
        [useless,unique_rows,useless2]=unique(round(start_from_genotype(:,1:end-2)*100),'rows');
        start_from_genotype=start_from_genotype(unique_rows,:);
        [useless,sorted_rows]=sort(start_from_genotype(:,end-1),'descend');
        start_from_genotype=start_from_genotype(sorted_rows,:);
        effective_childen_num=size(start_from_genotype,1);
        disp(['Effective population size=',num2str(effective_childen_num)])
        if effective_childen_num<population_size
            childen_genotype_temp=start_from_genotype;
        else
            childen_genotype_temp=start_from_genotype(1:population_size,:);
        end
        childen_genotype_list=childen_genotype_temp(:,1:end-2);
        sorted_index=1:size(childen_genotype_list,1);
        GW_sorted=childen_genotype_temp(:,end-1)';
        GN_sorted=childen_genotype_temp(:,end)';
    end
end

fnout=fopen('genotype_evolution_GW_GN.txt','a');
fnout2=fopen('genotype_best.txt','a');

if exist('genotype_evolution_GW_GN.txt', 'file')~=2 || exist_childen_num<population_size
    c_list=zeros(0,49);
    fnout_list={};
    fnout_list{population_size-exist_childen_num}='temp';
    for i=exist_childen_num+1:population_size
        individual_geno=zeros(1,47);
        for j=1:47
            coef_temp=rand();
            individual_geno(j)=PARA_list(j,1)*coef_temp + PARA_list(j,2)*(1-coef_temp);
            individual_geno(j)=round(individual_geno(j)*100)/100;
        end
        stem_NSC_Cost=((1-individual_geno(30))*Root_DW + (1-individual_geno(17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + (1-individual_geno(3))*Leaf_DW) ...
                            * 1.3     *     (1 - 0.03) / stem_starch; % relative stem NSC gain/loss
        if 1+stem_NSC_Cost<0 % restricting root/leaf/ear gain to keep stem NSC>=0
            scale_coef=-1/stem_NSC_Cost;
            individual_geno(30)=1-scale_coef+individual_geno(30)*scale_coef;
            individual_geno(17)=1-scale_coef+individual_geno(17)*scale_coef;
            individual_geno(3)=1-scale_coef+individual_geno(3)*scale_coef;
			stem_NSC_Cost=((1-individual_geno(30))*Root_DW + (1-individual_geno(17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + (1-individual_geno(3))*Leaf_DW) ...
							* 1.3     *     (1 - 0.03) / stem_starch;
        elseif 1+stem_NSC_Cost>stem_maxVolume/stem_starch % restricting root/leaf/ear loss to keep stem NSC<=stem_maxVolume
            scale_coef=(stem_maxVolume/stem_starch-1)/stem_NSC_Cost;
            individual_geno(30)=1-scale_coef+individual_geno(30)*scale_coef;
            individual_geno(17)=1-scale_coef+individual_geno(17)*scale_coef;
            individual_geno(3)=1-scale_coef+individual_geno(3)*scale_coef;
			stem_NSC_Cost=((1-individual_geno(30))*Root_DW + (1-individual_geno(17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + (1-individual_geno(3))*Leaf_DW) ...
							* 1.3     *     (1 - 0.03) / stem_starch;
        end

        delta_Area_leaf=(individual_geno(3)-1)*Area_leaf;
		delta_Vol_leaf=delta_Area_leaf*Thick_leaf*Vol_Fresh_ratio_leaf;
		delta_ON_IN_leaf=(conc_organic_N_leaf+conc_inorganic_N_leaf)*delta_Vol_leaf*14*6.25;
		delta_struN_leaf=delta_Area_leaf/SLA0*leaf_structural_N*6.25;
		conc_leaf_pro_new=(conc_leaf_pro*Leaf_DW-delta_ON_IN_leaf-delta_struN_leaf)/(Leaf_DW+delta_Area_leaf/SLA0);
		individual_geno(44)=conc_leaf_pro_new/conc_leaf_pro;
		delta_struN_grain=(individual_geno(17)-1)*Grain_NUM*Grain_VOL/1000.0/0.8*0.2*leaf_structural_N*6.25;
		delta_root_DW=(individual_geno(30)-1)*Root_DW;
		delta_struN_root=delta_root_DW*root_structural_N*6.25;
		delta_Vol_root=delta_root_DW/Fresh_dry_ratio_root*Vol_Fresh_ratio_root*10^-6;
		delta_ON_IN_root=(conc_organic_N_root+conc_inorganic_N_root)*delta_Vol_root*14*6.25;
		delta_stem_pro=-delta_struN_grain-delta_struN_root-delta_ON_IN_root;
		stem_N_Cost=delta_stem_pro/Stem_pro;
        if 1+stem_N_Cost<0 % restricting root/leaf/ear gain to keep stem N>=0
            scale_coef=-1/stem_N_Cost;
            individual_geno(30)=1-scale_coef+individual_geno(30)*scale_coef;
            individual_geno(17)=1-scale_coef+individual_geno(17)*scale_coef;
			delta_struN_grain=(individual_geno(17)-1)*Grain_NUM*Grain_VOL/1000.0/0.8*0.2*leaf_structural_N*6.25;
			delta_root_DW=(individual_geno(30)-1)*Root_DW;
			delta_struN_root=delta_root_DW*root_structural_N*6.25;
			delta_Vol_root=delta_root_DW/Fresh_dry_ratio_root*Vol_Fresh_ratio_root*10^-6;
			delta_ON_IN_root=(conc_organic_N_root+conc_inorganic_N_root)*delta_Vol_root*14*6.25;
			delta_stem_pro=-delta_struN_grain-delta_struN_root-delta_ON_IN_root;
			stem_NSC_Cost=((1-individual_geno(30))*Root_DW + (1-individual_geno(17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + (1-individual_geno(3))*Leaf_DW) ...
							* 1.3     *     (1 - 0.03) / stem_starch;
			stem_N_Cost=delta_stem_pro/Stem_pro;
        end

		individual_geno(47)=max(0,1+stem_N_Cost);
        individual_geno(25)=max(0,1+stem_NSC_Cost);

        c_list(i-exist_childen_num,:)=[[i,i],individual_geno];
        fnout_list{i-exist_childen_num}=strcat(num2str(i),'-plot.txt');
    end
    
    fnCommand_temp=fopen('batchCommands.sh','w');
    for i=1:length(fnout_list) %% multi-threading
        screen_info=['echo the ',num2str(round(c_list(i,1))),'th running...\n'];
        fprintf(fnCommand_temp,screen_info);
        inputString=num2str(c_list(i,:),'%.6f,');
        inputString=['[',inputString(1:end-1),']'];
        command_temp=['matlab -nodesktop -singleCompThread -nojvm -nosplash -nodisplay -r "RunModel(',inputString,',0,0,0,0);exit" &\n'];
        fprintf(fnCommand_temp,command_temp);
        if ~mod(i,maxPoolNum)
            fprintf(fnCommand_temp,'wait\n');
        end
    end
    fprintf(fnCommand_temp,'wait\n');
    fclose(fnCommand_temp);
    disp('Batch run begins ...')
    system('sh batchCommands.sh');
    
    childen_num=size(c_list,1);
    Total_Biomass_solid_soluble=zeros(1,childen_num);
    N_grain_content=zeros(1,childen_num);
    for i=exist_childen_num+1:population_size
        try
            yied_data=importdata(fnout_list{i-exist_childen_num});
            yied_data=yied_data.data;
            Total_Biomass_solid_soluble(i-exist_childen_num)=yied_data(end,22);
            N_grain_content(i-exist_childen_num)=yied_data(end,25);
            content=[c_list(i-exist_childen_num,:),Total_Biomass_solid_soluble(i-exist_childen_num),N_grain_content(i-exist_childen_num)];
            fprintf(fnout,'%.2f\t',content);
            fprintf(fnout,'\n');
        catch
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            disp([fnout_list{i-exist_childen_num},': ','failed...']);
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        end
    end
    
    childen_genotype_list_temp=c_list;
    childen_genotype_list=[childen_genotype_list;childen_genotype_list_temp];    
    GW_sorted=[GW_sorted,Total_Biomass_solid_soluble];
    GN_sorted=[GN_sorted,N_grain_content];
    [GW_sorted,sorted_index]=sort(GW_sorted,'descend');
    GN_sorted=GN_sorted(sorted_index);
end

generation_now=max(1,floor(exist_childen_num/population_size));

%% evolution 
parents_genotype_WN_list=zeros(0,47); % gennotype,grain weight, grain nitrogen concentration
for generation_num=generation_now:generation
    try
        childen_genotype_list=[childen_genotype_list(sorted_index,:),GW_sorted',GN_sorted'];
    catch
        size(childen_genotype_list)
        size(sorted_index)
        size(GW_sorted)
        size(GN_sorted)
        error('ctg')
    end
    candidate_parents_genotype_WN_list=[childen_genotype_list;parents_genotype_WN_list];
    [useless,unique_rows,useless2]=unique(round(candidate_parents_genotype_WN_list(:,3:end-2)*1000),'rows');
    candidate_parents_genotype_WN_list=candidate_parents_genotype_WN_list(unique_rows,:);
    [useless,index_WN]=sort(candidate_parents_genotype_WN_list(:,end-1),'descend'); % sorted by grain weight
    candidate_parents_genotype_WN_list=candidate_parents_genotype_WN_list(index_WN,:);
    parents_genotype_WN_list=candidate_parents_genotype_WN_list(1:parents_num,:); % update parents genotype
    GW_lowest=min(candidate_parents_genotype_WN_list(:,end-1));
    %GW_lowest=2;
    fit_sum = sum(candidate_parents_genotype_WN_list(:,end-1)-GW_lowest); % fitness function
    probability_list = (candidate_parents_genotype_WN_list(:,end-1)-GW_lowest)./fit_sum*300;
    probability_list(1:parents_num)=probability_list(1:parents_num)*3; % *3 to award the best individuals
    Russian_roulette_list=[];
    for temp_i = 1:length(unique_rows)
        Russian_roulette_list=[Russian_roulette_list,temp_i*ones(1,round(probability_list(temp_i)))];
    end
    Russian_roulette_size=length(Russian_roulette_list);
    fprintf(fnout2,'%.2f\t',childen_genotype_list(1,:));% recording best individual in current generation
    fprintf(fnout2,'\n');
    disp(['The ',num2str(generation_num),'th generation done! GW and GN of the best individual: ',num2str(childen_genotype_list(1,end-1)),', ',num2str(childen_genotype_list(1,end))])
    
    parents_genotype_list=candidate_parents_genotype_WN_list(:,1:end-2);
    childen_genotype_list=zeros(population_size+3,49);
    for i=1:(population_size+3)/4
        parents_temp=round(rand(1,2)*Russian_roulette_size);
        if parents_temp(1)*parents_temp(2)<0.1
            if parents_temp(1)<0.1
                parents_temp(1)=Russian_roulette_size;
            end
            if parents_temp(2)<0.1
                parents_temp(2)=Russian_roulette_size;
            end
        end
        parents_temp=Russian_roulette_list(parents_temp);
        recombination_pos=ceil(rand()*41);
        hybrid_children=zeros(4,47);
        
        hybrid_children(1,:)=[parents_genotype_list(parents_temp(1),3:3+recombination_pos),parents_genotype_list(parents_temp(2),3+recombination_pos+1:end)];
        hybrid_children(2,:)=[parents_genotype_list(parents_temp(2),3:3+recombination_pos),parents_genotype_list(parents_temp(1),3+recombination_pos+1:end)];
        hybrid_children(3,:)=(parents_genotype_list(parents_temp(1),3:end)+parents_genotype_list(parents_temp(2),3:end))/2;
        hybrid_children(3,:)=round(hybrid_children(3,:)*100)/100;
        hybrid_children(4,:)=max(parents_genotype_list(parents_temp(1),3:end),parents_genotype_list(parents_temp(2),3:end));
        for j=1:37*mutation_rate/100*5
            for child_k=1:4
                if rand()<0.2
                    mutation_pos=ceil(rand()*47);
                    hybrid_children(child_k,mutation_pos)=PARA_list(mutation_pos,1)+rand()*(PARA_list(mutation_pos,2)-PARA_list(mutation_pos,1));
                    hybrid_children(child_k,mutation_pos)=round(hybrid_children(child_k,mutation_pos)*100)/100;
                end
            end
        end
        for child_k=1:4
            stem_NSC_Cost=((1-hybrid_children(child_k,30))*Root_DW + (1-hybrid_children(child_k,17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + ...
                (1-hybrid_children(child_k,3))*Leaf_DW) * 1.3     *     (1 - 0.03) / stem_starch;
            if 1+stem_NSC_Cost<0
                scale_coef=-1/stem_NSC_Cost;
                hybrid_children(child_k,30)=1-scale_coef+hybrid_children(child_k,30)*scale_coef;
                hybrid_children(child_k,17)=1-scale_coef+hybrid_children(child_k,17)*scale_coef;
                hybrid_children(child_k,3)=1-scale_coef+hybrid_children(child_k,3)*scale_coef;
                stem_NSC_Cost=((1-hybrid_children(child_k,30))*Root_DW + (1-hybrid_children(child_k,17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + ...
                    (1-hybrid_children(child_k,3))*Leaf_DW) * 1.3     *     (1 - 0.03) / stem_starch;
            elseif 1+stem_NSC_Cost>stem_maxVolume/stem_starch
                scale_coef=(stem_maxVolume/stem_starch-1)/stem_NSC_Cost;
                hybrid_children(child_k,30)=1-scale_coef+hybrid_children(child_k,30)*scale_coef;
                hybrid_children(child_k,17)=1-scale_coef+hybrid_children(child_k,17)*scale_coef;
                hybrid_children(child_k,3)=1-scale_coef+hybrid_children(child_k,3)*scale_coef;
                stem_NSC_Cost=((1-hybrid_children(child_k,30))*Root_DW + (1-hybrid_children(child_k,17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + ...
                    (1-hybrid_children(child_k,3))*Leaf_DW) * 1.3     *     (1 - 0.03) / stem_starch;
            end
            delta_Area_leaf=(hybrid_children(child_k,3)-1)*Area_leaf;
            delta_Vol_leaf=delta_Area_leaf*Thick_leaf*Vol_Fresh_ratio_leaf;
            delta_ON_IN_leaf=(conc_organic_N_leaf+conc_inorganic_N_leaf)*delta_Vol_leaf*14*6.25;
            delta_struN_leaf=delta_Area_leaf/SLA0*leaf_structural_N*6.25;
            conc_leaf_pro_new=(conc_leaf_pro*Leaf_DW-delta_ON_IN_leaf-delta_struN_leaf)/(Leaf_DW+delta_Area_leaf/SLA0);
            hybrid_children(child_k,44)=conc_leaf_pro_new/conc_leaf_pro;
            delta_struN_grain=(hybrid_children(child_k,17)-1)*Grain_NUM*Grain_VOL/1000.0/0.8*0.2*leaf_structural_N*6.25;
            delta_root_DW=(hybrid_children(child_k,30)-1)*Root_DW;
            delta_struN_root=delta_root_DW*root_structural_N*6.25;
            delta_Vol_root=delta_root_DW/Fresh_dry_ratio_root*Vol_Fresh_ratio_root*10^-6;
            delta_ON_IN_root=(conc_organic_N_root+conc_inorganic_N_root)*delta_Vol_root*14*6.25;
            delta_stem_pro=-delta_struN_grain-delta_struN_root-delta_ON_IN_root;
            stem_N_Cost=delta_stem_pro/Stem_pro;
            if 1+stem_N_Cost<0
                scale_coef=-1/stem_N_Cost;
                hybrid_children(child_k,30)=1-scale_coef+hybrid_children(child_k,30)*scale_coef;
                hybrid_children(child_k,17)=1-scale_coef+hybrid_children(child_k,17)*scale_coef;
                delta_struN_grain=(hybrid_children(child_k,17)-1)*Grain_NUM*Grain_VOL/1000.0/0.8*0.2*leaf_structural_N*6.25;
                delta_root_DW=(hybrid_children(child_k,30)-1)*Root_DW;
                delta_struN_root=delta_root_DW*root_structural_N*6.25;
                delta_Vol_root=delta_root_DW/Fresh_dry_ratio_root*Vol_Fresh_ratio_root*10^-6;
                delta_ON_IN_root=(conc_organic_N_root+conc_inorganic_N_root)*delta_Vol_root*14*6.25;
                delta_stem_pro=-delta_struN_grain-delta_struN_root-delta_ON_IN_root;
                stem_NSC_Cost=((1-hybrid_children(child_k,30))*Root_DW + (1-hybrid_children(child_k,17))*Grain_NUM*Grain_VOL/1000.0/0.8*0.2 + ...
                    (1-hybrid_children(child_k,3))*Leaf_DW) * 1.3     *     (1 - 0.03) / stem_starch;
                stem_N_Cost=delta_stem_pro/Stem_pro;
            end
            hybrid_children(child_k,47)=max(0,1+stem_N_Cost);
            hybrid_children(child_k,25)=max(0,1+stem_NSC_Cost);
        end
        
        childen_genotype_list(4*i-3,:)=[population_size*generation_num+4*i-3,population_size*generation_num+4*i-3,hybrid_children(1,:)];
        childen_genotype_list(4*i-2,:)=[population_size*generation_num+4*i-2,population_size*generation_num+4*i-2,hybrid_children(2,:)];
        childen_genotype_list(4*i-1,:)=[population_size*generation_num+4*i-1,population_size*generation_num+4*i-1,hybrid_children(3,:)];
        childen_genotype_list(4*i,:)=[population_size*generation_num+4*i,population_size*generation_num+4*i,hybrid_children(4,:)];
    end
    childen_genotype_list=childen_genotype_list(1:population_size,:);
    fnout_list={};
    fnout_list{population_size}='temp';
    for i=population_size*generation_num+1:population_size*(generation_num+1)
        fnout_list{i-population_size*generation_num}=strcat(num2str(i),'-plot.txt');
    end
    
    childen_genotype_list_new=zeros(0,49);
    fnout_list_new={};
    for i=1:population_size
        if exist(fnout_list{i},'file')~=2 % don't rerun for existing files
            childen_genotype_list_new=[childen_genotype_list_new;childen_genotype_list(i,:)];
            fnout_list_new{end+1}=fnout_list{i};
        else
            disp(['the ',num2str(population_size*generation_num+i),'th running...']);
        end
    end
    if isempty(fnout_list_new)
        continue
    end
    
    childen_genotype_list=childen_genotype_list_new;
    fnout_list=fnout_list_new;
    fnCommand_temp=fopen('batchCommands.sh','w');
    for i=1:length(fnout_list) %% multi-threading
        screen_info=['echo the ',num2str(round(childen_genotype_list(i,1))),'th running...\n'];
        fprintf(fnCommand_temp,screen_info);
        inputString=num2str(childen_genotype_list(i,:),'%.6f,');
        inputString=['[',inputString(1:end-1),']'];
        command_temp=['matlab -nodesktop -singleCompThread -nojvm -nosplash -nodisplay -r "RunModel(',inputString,',0,0,0,0);exit" &\n'];
        fprintf(fnCommand_temp,command_temp);
        if ~mod(i,maxPoolNum)
            fprintf(fnCommand_temp,'wait\n');
        end
    end
    fprintf(fnCommand_temp,'wait\n');
    fclose(fnCommand_temp);
    disp('Batch run begins ...')
    system('sh batchCommands.sh');
    
    Total_Biomass_solid_soluble=[];
    N_grain_content=[];
    for i=1:length(fnout_list)
        try
            yied_data=importdata(fnout_list{i});
            yied_data=yied_data.data;
            Total_Biomass_solid_soluble(end+1)=yied_data(end,22);
            N_grain_content(end+1)=yied_data(end,25);
            content=[childen_genotype_list(i,:),Total_Biomass_solid_soluble(end),N_grain_content(end)];
            fprintf(fnout,'%.2f\t',content);
            fprintf(fnout,'\n');
        catch
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            disp([fnout_list{i},': ','failed...']);
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        end
    end
    
    [GW_sorted,sorted_index]=sort(Total_Biomass_solid_soluble,'descend');
    GN_sorted=N_grain_content(sorted_index);
end

%% write the last generation info. into file
Total_Biomass_solid_soluble=zeros(1,population_size);
N_grain_content=zeros(1,population_size);
for i=1:population_size
    try
        yied_data=importdata(fnout_list{i});
        yied_data=yied_data.data;
        Total_Biomass_solid_soluble(i)=yied_data(end,22);
        N_grain_content(i)=yied_data(end,25);
    catch
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        disp([fnout_list{i},': ','failed...']);
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    end
end
[GW_sorted,sorted_index]=sort(Total_Biomass_solid_soluble,'descend');
GN_sorted=N_grain_content(sorted_index);
for i=1:population_size
    content=[childen_genotype_list(i,:),Total_Biomass_solid_soluble(i),N_grain_content(i)];
    fprintf(fnout,'%.2f\t',content);
    fprintf(fnout,'\n');
end
childen_genotype_list=[childen_genotype_list(sorted_index,:),GW_sorted',GN_sorted'];
fprintf(fnout2,'%.2f\t',childen_genotype_list(1,:));
fprintf(fnout2,'\n');
disp(['The ',num2str(generation),'th generation done! GW and GN of the best individual: ',num2str(childen_genotype_list(1,end-1)),', ',num2str(childen_genotype_list(1,end))])
fclose(fnout);
fclose(fnout2);
