function [Tt,simulation]=RunModel(c,shading_on,shade_period,plant_spacing_on,CN_economy_print_on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RunModel.m: This function defines plant initial state at flowering, calls simulation_main.m  to simulate plant physiolocial change from flowering to harvest,
%             and summarizes results into text file.
%
%  Input:
%   c: A numerical array or cell struc with 49 elements for setting adjustment coefficients of model parameters
%   shading_on: If shading_on=1, light intensity during the first shade_period days after flowering will be halved
%   shade_period: The time period for shading (days)
%   plant_spacing_on: If plant_spacing_on=1, planting density after shading will be halved
%   CN_economy_print_on: If CN_economy_print_on=1, plant carbon and nitrogen mass balance and economy during grain filling will be printed on the screen
%
%  Output:
%   Tt: Time points from flowering to harvest (days)
%   simulation: Raw simulation result matrix (values of traits at each Tt)
%   c{1}.txt: raw results regarding to plant physiolocial change during grain filling
%   c{1}_plot.txt: processed and selected results used for ploting
%   If CN_economy_print_on=1, plant carbon and nitrogen mass balance and economy during grain filling will be printed on the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic
    global shading_span;
    shading_span=shade_period;
    global SPACING;
    SPACING=plant_spacing_on;
    main_dir='./';
    %% pertubation study
    if iscell(c)
        fnOut_basic=[num2str(c{1}),'.txt']; % simulation results, single tiller based
        fnOut_plot=[num2str(c{2}),'-plot.txt']; % selected/processed simulation results, single tiller based
        c=cell2mat(c(3:end));
        c=[1,1,c];
    else
        fnOut_basic=[num2str(c(1)),'.txt'];
        fnOut_plot=[num2str(c(2)),'-plot.txt'];
    end
    %% input coef. of variables 
    coef_Area_leaf=c(5);
    coef_grainNum=c(19);
    coef_grain_wid_max=c(25);
    coef_grain_len_wid_ratio=c(26);
    coef_Star_stem=c(27);
    coef_Fresh_root=c(32);
    coef_conc_iN_soil=c(33);
    coef_R_l_s_phloem0=c(41);
    coef_R_r_l_phloem0=c(42);
    coef_tillerNum=c(44);
    coef_CO2a=c(45);
    coef_Pro_content_leaf=c(46);
    coef_grain_struc_N=c(47);
    coef_stem_struc=c(48);
    coef_Pro_stem=c(49);
    
    %% constants
    start_time=0;
    end_time=86400*50; % Unit: seconds
    grow_efficiency=0.77; 
    grow_efficiency_grain=0.77; 
    tillerNum=coef_tillerNum*12; % tillers per plant
    plantDens=25; % plants per m2
    coef_Suc_maintain_leaf=6.9*10^-5;
    coef_Suc_maintain_grain=6.9*10^-5;
    coef_Suc_maintain_root=1.39*10^-4;
    coef_Suc_maintain_stem=6.9*10^-5;
    Suc_molWeight=342;
    N_molWeight=14;
    Pro_molWeight=87.5; % 14*6.25
    N_to_Pro_mass=6.25;
    
    %% global input variables from files
    weather_data=importdata([main_dir,'/weather_input.txt']);
    global PPFD_list;
    PPFD_list=weather_data.data(:,2); % umol/m2/s
    if shade_period==1
        PPFD_list(1:240,1)=PPFD_list(1:240,1)*(100-shading_on)/100;
    elseif shade_period==2
        PPFD_list(241:480,1)=PPFD_list(241:480,1)*(100-shading_on)/100;
    end
    global VPD_leaf_list;
    VPD_leaf_list=weather_data.data(:,3); % kPa

    %% leaf (big-leaf model)
    Vol_Fresh_ratio_leaf=0.25;
    Thick_leaf=0.000175; % m 
    CO2_a=400*coef_CO2a; % ppm, ambient CO2 conc.
    lat_Shanghai=0.541; % 31/180*pi=0.541
    DOY_start=228; % 2016-8-15 is the 228th day of the year
    Area_leaf=coef_Area_leaf*200*10^-4; %m2
    Area_leaf0=Area_leaf;
    SLA0=0.02; % g/m2
    conc_leaf_Pro0=coef_Pro_content_leaf*0.035*N_to_Pro_mass; %  g/g 
    Vol_leaf=Area_leaf*Thick_leaf*Vol_Fresh_ratio_leaf; % m3
    Vol_leaf0=Vol_leaf;
    conc_Suc_leaf=50; % mol/m3
    conc_oN_leaf=100; % mol/m3
    conc_iN_leaf=50; % mol/m3
    totMol_Pro_leaf=Area_leaf0/SLA0*conc_leaf_Pro0/Pro_molWeight; % mol Pro
    totMol_Pro_leaf0=totMol_Pro_leaf; % mol Pro
    Star_leaf=3*10^-6; % mol Suc.
    conc_TP_leaf=4; % mol/m3

    %% root
    Vol_soil=0.001; % m3; 0.30 m * 1/(tillerNum*plantingDensity) m2 
    Rad_root=0.8*10^-4; %m 
    Vol_Fresh_ratio_root=0.25;
    Fresh_dry_ratio_root=0.1; % convert root fresh weight to dry weight
    Absorb_area_ratio=1;
    Soil_water_content=0.25; %
    Vol_soil_water=Vol_soil*Soil_water_content; % m3
    conc_iN_soil=coef_conc_iN_soil*0.14; % mol/m3, equal to 2mg/L
    Fresh_root=coef_Fresh_root*5;
    Vol_root=Vol_Fresh_ratio_root*Fresh_root*10^-6; % m3
    Dry_root=Fresh_dry_ratio_root*Fresh_root;
    Vol_root0=Vol_root;
    Area_root=1/Vol_Fresh_ratio_root*2*Vol_root/Rad_root; % m2
    Area_root0=Area_root;
    conc_Suc_root=20; % mol/m3
    conc_iN_root=20; % mol/m3
    conc_oN_root=10; % mol/m3
    root_struc_N=1/200; %

    %% stem
    Vol_xylem=0.1*10^-6; % m3
    Vol_rphloem=0.03*10^-6; % m3
    Vol_lphloem=0.03*10^-6; % m3
    Vol_sphloem=0.03*10^-6; % m3
    R_l_s_phloem0=coef_R_l_s_phloem0*1.75*10^13;
    R_r_l_phloem0=coef_R_r_l_phloem0*7*10^13;
    conc_iN_xylem=40; % mol/m3
    conc_oN_xylem=10; % mol/m3
    conc_oN_rphloem=150; % mol/m3
    conc_Suc_rphloem=300; % mol/m3
    conc_oN_lphloem=200; % mol/m3
    conc_Suc_lphloem=600; % mol/m3
    conc_oN_sphloem=150; % mol/m3
    conc_Suc_sphloem=500; % mol/m3

    %% ear
    Vol_Fresh_ratio_grain=0.9; % Converting fresh weight to grain volume
    Fresh_new_dry_ratio_grain=0.1; % Converting fresh weight to dry weight
    pho_Pro=1; %t/m3, density of Pro
    pho_Star=1; %t/m3, density of Pro
    GSR=0.875; % grain setting rate
    grainNum=coef_grainNum*160*GSR;
    grain_wid_max=coef_grain_wid_max*2.4*10^-3;
    grain_len_wid_ratio=coef_grain_len_wid_ratio*2.8;
    Vol_grain_max=2.0267e-08*coef_grain_wid_max^3*coef_grain_len_wid_ratio; %m3/grain, max volume of single grain
    Area_all_grain_max=grainNum*4*pi*((1+2*grain_len_wid_ratio^1.6075)/3)^(1/1.6075)*(grain_wid_max/2)^2; % m2/ear, ellipsoid model
    Vol_grain=Vol_grain_max*0.001; % m3, initial volume of a single grain
    Area_grain=4*pi*((1+2*grain_len_wid_ratio^1.6075)/3)^(1/1.6075) ...
        *(Vol_grain*3/(4*pi)/grain_len_wid_ratio)^(2/3); % m2 (single grain), ellipsoid model, S=4pi*((1+2*L^1.6075)/3)^(1/1.6075)*r2, V=4/3*pi*r3*L
    conc_Suc_grain=250; % mol/m3, initial grain Suc concentration
    conc_oN_grain=50; % mol/m3
    store_Pro_ear=0; % g DW/ear
    store_Star_ear=0; % g DW/ear
    grain_struc_N=coef_grain_struc_N*4/100;
    
    %% stem
    FW_DW_stem=3; % Converting dry weight to stem fresh weight 
    W_stem_struc=coef_stem_struc*3; % g
    Vol_stem=W_stem_struc*FW_DW_stem*0.2*10^-6; %m3, assume 20% of the stem volume is stem
    conc_Suc_stem=conc_Suc_lphloem;
    conc_oN_stem=conc_oN_lphloem;
    Pro_stem=coef_Pro_stem*0.02; % g DW
    Star_stem=coef_Star_stem*0.75; % g DW
    Vol_soluble_stem=(Vol_stem-Pro_stem*10^-6/pho_Pro-Star_stem*10^-6/pho_Star); % m3
    
    %% record variable
    % root
    Suc_in_root=0;
    oN_shoot_root_netFlux=0;
    iN_out_root=0;
    Vol_sene_root=0;
    Vol_grow_root=0;
    Rd_root=0;
    N_assimilation_root=0;
    N_absorb_root=0;
    N_shoot_root_netFlux=0;
    oN_import_root=0;
    oN_export_root=0;
    iN_export_root=0;
    % leaf
    CO2_assimilation_leaf=0;
    Suc_out_leaf=0;
    oN_out_leaf=0;
    iN_in_leaf=0;
    oN_in_leaf=0;
    N_in_leaf=0;
    N_ass_leaf=0;
    N_leak_leaf=0;
    C_leak_leaf=0;
    Area_sene_leaf=0;
    Rd_leaf=0;
    gs_leaf=0;
    CO2i=0;
    % ear
    Suc_in_ear=0;
    oN_in_ear=0;
    Vol_grow_ear=0;
    Biomass_grow_ear=0;
    Rd_ear=0;
    % flow
    Flow_lphloem_rphloem=0;
    Flow_lphloem_sphloem=0;
    oN_xylem_to_lphloem=0;
    % stem
    Rd_stem=0;
    
    %% packing for constants
    constants_list=[Suc_molWeight,N_molWeight,Pro_molWeight,N_to_Pro_mass,root_struc_N,grain_struc_N,grow_efficiency,...
        grow_efficiency_grain,tillerNum,plantDens,coef_Suc_maintain_leaf,coef_Suc_maintain_grain,coef_Suc_maintain_root,coef_Suc_maintain_stem,...%14
        Vol_Fresh_ratio_leaf,Thick_leaf,CO2_a,lat_Shanghai,DOY_start,Area_leaf0,Vol_leaf0,totMol_Pro_leaf0,...%8
        Vol_soil,Rad_root,Vol_Fresh_ratio_root,Fresh_dry_ratio_root,Absorb_area_ratio,Soil_water_content,Vol_root0,Area_root0,... %8
        Vol_xylem,Vol_rphloem,Vol_lphloem,Vol_sphloem,R_l_s_phloem0,R_r_l_phloem0,... %6
        Vol_Fresh_ratio_grain,Fresh_new_dry_ratio_grain,pho_Pro,pho_Star,grainNum,grain_len_wid_ratio,Vol_grain_max,Area_all_grain_max,... %8
        Vol_stem]; %1 , in total: 45
    root_input=[Vol_root,Dry_root,conc_Suc_root,conc_iN_root,conc_oN_root,conc_iN_soil,Vol_soil_water,zeros(1,100-45-7),...
        constants_list]'; 
    leaf_input=[Area_leaf,Vol_leaf,conc_Suc_leaf,conc_iN_leaf,conc_oN_leaf,totMol_Pro_leaf,Star_leaf,conc_TP_leaf,zeros(1,100-8)]'; 
    ear_input=[Vol_grain,Area_grain,conc_Suc_grain,conc_oN_grain,store_Pro_ear,store_Star_ear,zeros(1,100-6)]'; 
    flow_input=[conc_iN_xylem,conc_Suc_rphloem,conc_Suc_lphloem,conc_Suc_sphloem,conc_oN_rphloem,conc_oN_lphloem,...
        conc_oN_sphloem,conc_oN_xylem,zeros(1,100-8)]'; 
    stem_input=[conc_Suc_stem,conc_oN_stem,Star_stem,Pro_stem,Vol_soluble_stem,zeros(1,100-43-5),c(3:45)]'; 
    record_input=[Suc_in_root,oN_shoot_root_netFlux,iN_out_root,Vol_sene_root,Vol_grow_root,Rd_root,N_assimilation_root,N_absorb_root,...
        N_shoot_root_netFlux,oN_import_root,oN_export_root,iN_export_root,zeros(1,20-12),... % root:12
        CO2_assimilation_leaf,Suc_out_leaf,oN_out_leaf,iN_in_leaf,oN_in_leaf,N_in_leaf,N_ass_leaf,N_leak_leaf,C_leak_leaf,...
        Area_sene_leaf,Rd_leaf,gs_leaf,CO2i,zeros(1,20-13),... % leaf:13
        Suc_in_ear,oN_in_ear,Vol_grow_ear,Biomass_grow_ear,Rd_ear,zeros(1,20-5),... % grain:5
        Flow_lphloem_rphloem,Flow_lphloem_sphloem,oN_xylem_to_lphloem,zeros(1,20-3),... % flow:3
        Rd_stem,zeros(1,20-1)]'; % stem:1
    options= odeset('Events',@MyEventFunction);
    try
        [Tt, simulation]=ode15s(@simulation_main, [start_time:3600:end_time], [root_input,leaf_input,ear_input,flow_input,stem_input,record_input],options);
    catch
        disp('Time Out!');
        Tt=[];
        simulation=[];
        return
    end
    % convert results from tiller to plant
    columns=[1,2,7,...
             101,102,106,107,...
             205,206,...
             403:405,...
             501:512,521:531,541:545,561:563,581];
    simulation(:,columns)=simulation(:,columns)*tillerNum;
    
    %% write to txt (*.txt and *-plot.txt)
    %%%%%%%*.txt%%%%%%%
    fnout=fopen([main_dir,fnOut_basic],'w');
    fprintf(fnout,[...
        'Time (s)\tVolume_root (m3)\tWeight_root (g)\t[Suc_root] (mol/m3)\t[iN_root] (mol/m3)\t[oN_root] (mol/m3)\t[iN_soil] (mol/m3)\t'...
        'Volume_soil_water (m3)\t'... % root: 7
        'Area_leaf (m2)\tVolume_leaf (m3)\t[Suc_leaf] (mol/m3)\t[iN_leaf] (mol/m3)\t[oN_leaf] (mol/m3)\tPro_leaf (mol)\t'...
        'Star_leaf (mol Equiv. Suc)\t[TP_leaf] (mol Equiv. Suc/m3)\t'... % leaf: 8
        'Volume_grain (m3)\tArea_grain (m2)\t[Suc_grain] (mol/m3)\t[oN_grain] (mol/m3)\tPro_ear (g)\tStar_ear (g)\t'... % ear: 6
        '[iN_xylem] (mol/m3)\t[Suc_rphloem] (mol/m3)\t[Suc_lphloem] (mol/m3)\t[Suc_sphloem] (mol/m3)\t[oN_rphloem] (mol/m3)\t'...
        '[oN_lphloem] (mol/m3)\t[oN_sphloem] (mol/m3)\t[oN_xylem] (mol/m3)\t'... % flow: 8
        '[Suc_stem] (mol/m3)\t[oN_stem] (mol/m3)\tStar_stem (g)\tPro_stem (g)\tVolume_soluble_stem (m3)\t'... % stem: 5
        'Suc_root_imported (mol)\toN_root_imported (mol)\tiN_root_exported (mol)\tVolume_root_senesced (m3)\tVolume_root_grew (m3)\t'...
        'Rd_root (mol Equiv. Suc)\tN_root_assimid (mol)\tN_root_absorbed (mol)\ttotalN_shoot2root_netFlux (mol)\toN_root_imported (mol)\t'...
        'oN_root_exported (mol)\tiN_root_exported (mol)\t'... % record_root: 12
        'CO2_canopy_assimid (mol Equiv. Suc)\tSuc_leaf_exported (mol)\toN_leaf_exported (mol)\tiN_leaf_imported (mol)\toN_leaf_imported (mol)\t'...
        'totalN_leaf_imported (mol)\tiN_leaf_assimi (mol)\tN_leaf_leaked (mol)\tC_leaf_leaked (mol Equiv. Suc)\tArea_leaf_senesced (m2)\t'...
        'Rd_leaf (mol Equiv. Suc)\t'... % record_leaf: 11
        'Suc_ear_imported (mol)\toN_ear_imported (mol)\tVolume_ear_grew (m3)\tBiomass_ear_grew (g)\tRd_ear (mol Equiv. Suc)\t'... % record_ear: 5
        'Fluxes_lphloem_rphloem (m3)\tFluxes_lphloem_sphloem (m3)\toN_xylem_to_lphloem (mol)\t'... % record_flow: 3
        'Rd_stem (mol Equiv. Suc)\n'... % record_stem: 1
                  ]);
    content=[Tt,simulation(:,[1:7,101:108,201:206,301:308,401:405,501:512,521:531,541:545,561:563,581])]'; % 67 items
    format=[repmat('%.2e\t',1,66),'%.2e\n'];
    fprintf(fnout,format,content);
    fclose(fnout);
    
    %%%%%%%*-plot.txt%%%%%%%
    timeDay=Tt/86400;
    % root
    Vol_root_Total=10^6*simulation(:,1); %cm3
    Vol_root_grow=10^6*simulation(:,505); %cm3
    Vol_root_sene=10^6*simulation(:,504); %cm3
    Rd_root=86400*10^6*diff(simulation(:,506))./diff(Tt)*12; % umol CO2/d
    Rd_root=[Rd_root(1);Rd_root];
    Unit_Rd_root=86400*10^6*diff(simulation(:,506))./diff(Tt)./(10^6*simulation(2:end,1)/Vol_Fresh_ratio_root)*12; % umol CO2/FW/d
    Unit_Rd_root=[Unit_Rd_root(1);Unit_Rd_root];
    N_absorb_root=10^6*simulation(:,508); % umol
    N_assimi_root=10^6*simulation(:,507); % umol
    conc_Suc_root=simulation(:,3); % mol/m3
    conc_oN_root=simulation(:,5); % mol/m3
    conc_iN_root=simulation(:,4); % mol/m3
    % leaf
    Area_leaf=10^4*simulation(:,101); % cm2
    Photosyn_canopy=86400*10^6*diff(simulation(:,521))./diff(Tt)*12; %umol CO2/d
    Photosyn_canopy=[Photosyn_canopy(1);Photosyn_canopy];
    Rd_leaf=86400*10^6*diff(simulation(:,531))./diff(Tt)*12; % umol CO2/d
    Rd_leaf=[Rd_leaf(1);Rd_leaf];
    NSC_leaf=(simulation(:,107)+simulation(:,103).*simulation(:,102))./simulation(:,101); % Star + Suc, mol Suc/m2
    conc_oN_leaf=simulation(:,105); % mol/m3
    conc_iN_leaf=simulation(:,104); % mol/m3
    pro_content_leaf=simulation(:,106); % mol
    % ear
    Vol_solid_single_grain=(simulation(:,205)/grainNum/tillerNum/pho_Pro+simulation(:,206)/grainNum/tillerNum/pho_Star)*10^-6; % m3/grain
    Struc_Biomass=simulation(:,544)-simulation(:,206)-simulation(:,205); % g/plant
    Total_Biomass_solid_soluble=simulation(:,544)+(simulation(:,201)-Vol_solid_single_grain)*grainNum*tillerNum.*(simulation(:,203)*Suc_molWeight ...
        + simulation(:,204)*Pro_molWeight); % g/plant
    Total_N_free_combined=simulation(:,205)/N_to_Pro_mass+Struc_Biomass*grain_struc_N+(simulation(:,201)-Vol_solid_single_grain)...
        * grainNum*tillerNum.*simulation(:,204)*N_molWeight; % g/plant
    Vol_ear=10^6*grainNum*tillerNum*simulation(:,201); % cm3/plant
    Vol_ear_soluble=10^6*grainNum*tillerNum*(simulation(:,201)-Vol_solid_single_grain); % cm3/plant
    Rd_ear=86400*10^6*diff(simulation(:,545))./diff(Tt)*12; %umol CO2/d
    Rd_ear=[Rd_ear(1);Rd_ear];
    Suc_grain=simulation(:,203).*(simulation(:,201)-Vol_solid_single_grain)*grainNum*tillerNum*Suc_molWeight*1000./Total_Biomass_solid_soluble; % mg/g
    oN_grain=simulation(:,204).*(simulation(:,201)-Vol_solid_single_grain)*grainNum*tillerNum*10^3./Total_Biomass_solid_soluble; % mmol/g
    N_grain_content=100*Total_N_free_combined./Total_Biomass_solid_soluble; % %
    % stem
    conc_Suc_rphloem=simulation(:,302); % mol/m3
    conc_Suc_lphloem=simulation(:,303); % mol/m3
    conc_Suc_sphloem=simulation(:,304); % mol/m3
    conc_oN_rphloem=simulation(:,305); % mol/m3
    conc_oN_lphloem=simulation(:,306); % mol/m3
    conc_oN_sphloem=simulation(:,307); % mol/m3
    conc_iN_xylem=simulation(:,301); % mol/m3
    conc_oN_xylem=simulation(:,308); % mol/m3
    conc_Suc_stem=simulation(:,401); % mol/m3
    Star_stem=simulation(:,403); % g
    Pro_stem=simulation(:,404); % g Pro
    T_Suc_in_root=simulation(:,501)*Suc_molWeight; % g
    Rd_stem=86400*10^6*diff(simulation(:,581))./diff(Tt)*12; %umol CO2/d
    Rd_stem=[Rd_stem(1);Rd_stem];
    
    fnout=fopen([main_dir,fnOut_plot],'w');
    fprintf(fnout,[...
        'timeDay (d)\tVolume_root (cm3)\tVolume_root_grew (cm3)\tVolume_root_senesced (cm3)\tRd_root (umol CO2/d)\t'...
        'Unit_Rd_root (umol CO2/g FreshWeight/s)\tN_absorb_root (umol)\tN_assimi_root (umol)\toN_root (mol/m3)\tSuc_root (mol/m3)\tiN_root (mol/m3)\t'... % root
        'Area_leaf (cm2)\tPhotosyn_canopy (umol Suc/d)\tRd_leaf (umol CO2/d)\tNSC_leaf (mol Suc/m2)\toN_leaf (mol/m3)\tiN_leaf (mol/m3)\t'...
        'pro_content_leaf (mol)\t'... % leaf
        'Volume_ear (cm3)\tVolume_ear_soluble (cm3)\tRd_ear (umol CO2/d)\tTotal_Biomass_solid_soluble (g/tiller)\tSuc_grain (mg/g)\toN_grain (umol/g)\t'...
        'N_grain_content (%%)\t'... % ear
        'Suc_rphloem (mol/m3)\tSuc_lphloem (mol/m3)\tSuc_sphloem (mol/m3)\toN_rphloem (mol/m3)\toN_lphloem (mol/m3)\toN_sphloem (mol/m3)\t'...
        'iN_xylem (mol/m3)\toN_xylem (mol/m3)\t'... % flow
        'Star_stem (g)\tT_Suc_in_root (g)\tconc_Suc_stem (mol/m3)\tRd_stem (umol CO2/d)\tPro_stem (g)\n'... % stem
       ]);
    data_out=[timeDay,Vol_root_Total,Vol_root_grow,Vol_root_sene,Rd_root,Unit_Rd_root,N_absorb_root,N_assimi_root,conc_oN_root,conc_Suc_root,conc_iN_root,... % root 11
              Area_leaf,Photosyn_canopy,Rd_leaf,NSC_leaf,conc_oN_leaf,conc_iN_leaf,pro_content_leaf,... % leaf 7
              Vol_ear,Vol_ear_soluble,Rd_ear,Total_Biomass_solid_soluble,Suc_grain,oN_grain,N_grain_content,... % ear 7
              conc_Suc_rphloem,conc_Suc_lphloem,conc_Suc_sphloem,conc_oN_rphloem,conc_oN_lphloem,conc_oN_sphloem,conc_iN_xylem,conc_oN_xylem,... % flow 8
              Star_stem,T_Suc_in_root,conc_Suc_stem,Rd_stem,Pro_stem]'; % stem 5. In total 38
    content=data_out;
    format=[repmat('%f\t',1,37),'%f\n'];
    fprintf(fnout,format,content);
    fclose(fnout);
    

    if CN_economy_print_on==1 % print plant carbon and nitrogen economy on screen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start: mass balance check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % root
        Suc_in_root=simulation(end,501); % mol
        N_abs_root=simulation(end,508);
        N_shoot_root_netFlux=simulation(end,509);
        Rd_root=simulation(end,506);
        NSC_root_flowering=simulation(1,1)*simulation(1,3); % mol Suc
        NSC_root_harvest=simulation(end,1)*simulation(end,3); % mol Suc
        N_root_flowering=simulation(1,1)*(simulation(1,4)+simulation(1,5)); % mol iN+oN
        root_grow_weight=simulation(end,505)*Fresh_dry_ratio_root/Vol_Fresh_ratio_root*10^6;
        N_root_harvest=simulation(end,1)*(simulation(end,4)+simulation(end,5)); % mol iN + oN
        N_root_harvest_struc= root_grow_weight*root_struc_N/N_molWeight;
        coef_C_skeleton_iN2oN=0.417;
        Suc_Nass_root=simulation(end,507)*coef_C_skeleton_iN2oN;
        Suc_grow_root=simulation(end,505)*Fresh_dry_ratio_root/Vol_Fresh_ratio_root*10^6*(1-root_struc_N*N_to_Pro_mass)/Suc_molWeight;
        delta_Suc_root=simulation(end,1)*simulation(end,3) - simulation(1,1)*simulation(1,3);
        delta_N_root=N_root_harvest - N_root_flowering;
        Root_C=[Suc_in_root-sum([Rd_root,Suc_Nass_root,Suc_grow_root])-delta_Suc_root,Suc_in_root,sum([Rd_root,Suc_Nass_root,Suc_grow_root]),...
            delta_Suc_root,Rd_root,Suc_Nass_root,Suc_grow_root]*1000; % mmol
        Root_N=[N_abs_root + N_shoot_root_netFlux - N_root_harvest_struc - delta_N_root, N_abs_root + N_shoot_root_netFlux, N_root_harvest_struc,...
            delta_N_root, simulation(end,1)*simulation(end,4), simulation(end,1)*simulation(end,5)]*1000; % mmol [iN, oN, Struc]
        
        % leaf
        Suc_in_leaf=-simulation(end,522);
        C_abs_leaf=simulation(end,521);
        N_in_leaf=simulation(end,526);
        xylem_oN_in_leaf=simulation(end,525);
        xylem_iN_in_leaf=simulation(end,524);
        phloem_oN_out_leaf=simulation(end,523);
        N_leak_leaf=simulation(end,528);
        C_leak_leaf=simulation(end,529);
        Rd_leaf=simulation(end,531);
        NSC_leaf_flowering=simulation(1,107)+(simulation(1,103)+simulation(1,108))*simulation(1,102); % mol Star + Suc + TP
        NSC_leaf_harvest=simulation(end,107)+(simulation(end,103)+simulation(end,108))*simulation(end,102); % mol Star + Suc + TP
                % 104: I-N; 105: O-N; 106: Pro; 
        N_leaf_flowering=simulation(1,106)+(simulation(1,104)+simulation(1,105))*simulation(1,102); % mol Pro + I-N + O-N
        N_leaf_harvest=simulation(end,106)+(simulation(end,104)+simulation(end,105))*simulation(end,102); % mol Pro + I-N + O-N
        delta_N_leaf=N_leaf_harvest - N_leaf_flowering;
        Suc_Nass_leaf=simulation(end,527)*coef_C_skeleton_iN2oN;
        delta_NSC_leaf=NSC_leaf_harvest-NSC_leaf_flowering;
        Leaf_C=[C_abs_leaf+Suc_in_leaf-Rd_leaf-Suc_Nass_leaf-C_leak_leaf-delta_NSC_leaf,C_abs_leaf,-Suc_in_leaf,...
            Rd_leaf+Suc_Nass_leaf+C_leak_leaf,delta_NSC_leaf,Rd_leaf,Suc_Nass_leaf,C_leak_leaf]*1000;
        Leaf_N=[N_in_leaf - delta_N_leaf - N_leak_leaf, N_in_leaf, N_leak_leaf, delta_N_leaf, xylem_iN_in_leaf,xylem_oN_in_leaf, phloem_oN_out_leaf,...
            simulation(end,106), simulation(end,104)*simulation(end,102), simulation(end,105)*simulation(end,102)]*1000; % mol [Pro, I-N, O-N]
             
        % ear
        Suc_in_ear=simulation(end,541);
        N_in_ear=simulation(end,542);
        Rd_ear=simulation(end,545);
        Star_ear_harvest=simulation(end,206)/Suc_molWeight;
        Pro_ear_harvest=simulation(end,205)/Pro_molWeight;
        Vol_solid_ear0=(simulation(1,205)/pho_Pro+simulation(1,206)/pho_Star)*10^-6;
        Suc_ear_flowering=(simulation(1,201)*grainNum*tillerNum - Vol_solid_ear0)*simulation(1,203);
        Vol_solid_ear1=(simulation(end,205)/pho_Pro+simulation(end,206)/pho_Star)*10^-6;
        Suc_ear_harvest=(simulation(end,201)*grainNum*tillerNum - Vol_solid_ear1)*simulation(end,203);
        struc_Suc_ear_harvest=(simulation(end,544)-simulation(end,206)-simulation(end,205))*(1-grain_struc_N*N_to_Pro_mass)/Suc_molWeight;
        NSC_ear_flowering=Suc_ear_flowering+simulation(1,206)/Suc_molWeight; % Suc + Star
        NSC_ear_harvest=Star_ear_harvest+Suc_ear_harvest+struc_Suc_ear_harvest;
        delta_NSC_ear=NSC_ear_harvest - NSC_ear_flowering;
        oN_ear_flowering=(simulation(1,201)*grainNum*tillerNum - Vol_solid_ear0)*simulation(1,204);
        oN_ear_harvest=(simulation(end,201)*grainNum*tillerNum - Vol_solid_ear1)*simulation(end,204);
        struc_oN_ear_harvest=(simulation(end,544)-simulation(end,206)-simulation(end,205))*grain_struc_N/N_molWeight;
        N_ear_flowering=oN_ear_flowering+simulation(1,205)/Pro_molWeight; % O-N + Pro
        N_ear_harvest=Pro_ear_harvest+oN_ear_harvest+struc_oN_ear_harvest;
        delta_N_ear=N_ear_harvest - N_ear_flowering;
        Ear_C=[Suc_in_ear - (delta_NSC_ear+Rd_ear),Suc_in_ear,sum([Rd_ear,Star_ear_harvest,struc_Suc_ear_harvest]),Suc_ear_harvest-Suc_ear_flowering,...
            Rd_ear,Star_ear_harvest,struc_Suc_ear_harvest]*1000;
        Ear_N=[N_in_ear - delta_N_ear, N_in_ear,struc_oN_ear_harvest,N_ear_harvest - struc_oN_ear_harvest - N_ear_flowering,...
            Pro_ear_harvest, oN_ear_harvest]*1000; % [Pro, oN, strucN]
        
        % stem
        Vol_xylem=Vol_xylem*tillerNum;
        Vol_rphloem=Vol_rphloem*tillerNum;
        Vol_lphloem=Vol_lphloem*tillerNum;
        Vol_sphloem=Vol_sphloem*tillerNum;
        Vol_soluble_stem0=simulation(1,405);
        Vol_soluble_stem1=simulation(end,405);
        stem_Pro0=Pro_stem(1)/Pro_molWeight;
        stem_Pro1=Pro_stem(end)/Pro_molWeight;
        stem_N0=simulation(1,402)*Vol_soluble_stem0+conc_oN_lphloem(1)*Vol_lphloem+conc_oN_rphloem(1)*Vol_rphloem+conc_oN_sphloem(1)*Vol_sphloem...
            + (conc_iN_xylem(1)+conc_oN_xylem(1))*Vol_xylem;
        stem_N1=simulation(end,402)*Vol_soluble_stem1+conc_oN_lphloem(end)*Vol_lphloem+conc_oN_rphloem(end)*Vol_rphloem+conc_oN_sphloem(end)*Vol_sphloem...
            + (conc_iN_xylem(end)+conc_oN_xylem(end))*Vol_xylem;
        N_stem_flowering=stem_N0+stem_Pro0;
        N_stem_harvest=stem_N1+stem_Pro1;
        delta_N_stem=N_stem_harvest - N_stem_flowering;
        N_in_stem  =-N_shoot_root_netFlux - N_in_ear - N_in_leaf;
        Vol_soluble_stem0=simulation(1,405);
        stem_Star0=Star_stem(1)/Suc_molWeight;
        stem_Suc0=simulation(1,401)*Vol_soluble_stem0 + conc_Suc_lphloem(1)*Vol_lphloem + conc_Suc_rphloem(1)*Vol_rphloem + conc_Suc_sphloem(1)*Vol_sphloem;
        NSC_stem_flowering=stem_Suc0+stem_Star0; % NSC in stem and phloem
        stem_Star1=Star_stem(end)/Suc_molWeight;
        Vol_soluble_stem1=simulation(end,405);
        stem_Suc1=simulation(end,401)*Vol_soluble_stem1 + conc_Suc_lphloem(end)*Vol_lphloem + conc_Suc_rphloem(end)*Vol_rphloem + conc_Suc_sphloem(end)*Vol_sphloem;
        NSC_stem_harvest=stem_Suc1+stem_Star1; % NSC in stem and phloem
        %NSC_stem_flowering=stem_Star0+simulation(1,401)*Vol_soluble_stem0;  % NSC in stem 
        %NSC_stem_harvest=stem_Star1+simulation(end,401)*Vol_soluble_stem1;  % NSC in stem 
        delta_NSC_stem=NSC_stem_harvest - NSC_stem_flowering;
        Suc_in_stem=-Suc_in_root-Suc_in_ear-Suc_in_leaf;
        Rd_stem=simulation(end,581);
        Stem_C=[Suc_in_stem - (delta_NSC_stem+Rd_stem),Suc_in_stem,Rd_stem,delta_NSC_stem,NSC_stem_harvest,NSC_stem_flowering]*1000; % mmol
        Stem_N=[delta_N_stem - N_in_stem, -N_shoot_root_netFlux - N_in_ear - N_in_leaf, 0,delta_N_stem, stem_N1, stem_Pro1]*1000; % [N, Pro]
        
        N_root_ass=simulation(end,507);
        oN_root_import=simulation(end,510);
        oN_root_export=simulation(end,511);
        iN_root_export=simulation(end,512);
        oN_leaf_export=simulation(end,523);
        iN_leaf_import=simulation(end,524);
        N_leaf_ass=simulation(end,527);
        oN_leaf_import=simulation(end,525);
        N_leaf_sene_leak=simulation(end,528);
        C_leaf_sene_leak=simulation(end,529);
        oN_xylem_phloem_transfer=simulation(end,563);
        other_N_budget=[N_root_ass,oN_root_import,oN_root_export,iN_root_export,oN_leaf_export,iN_leaf_import,N_leaf_ass,oN_leaf_import,N_leaf_sene_leak,...
            C_leaf_sene_leak,oN_xylem_phloem_transfer]*1000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end: mass balance check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('%%%%%%%%%%%%%%%%%%% start: mass balance check %%%%%%%%%%%%%%%%%%%')
        fprintf(['Root C balance (mmol/plant): [0, Suc import, Suc comsume, Suc surplus(Harvest-Heading), Rd, Nass, Growth]=[%.4f, %.4f, %.4f, %.4f, '...
            '%.4f, %.4f, %.4f]\n'],Root_C*12);
        fprintf(['Root N balance (mmol/plant): [0, N income(Absorb&Import), N consume(Growth), NSN surplus, Harvest(I-N), Harvest(O-N)]=[%.4f, %.4f, '...
            '%.4f, %.4f, %.4f, %.4f]\n'],Root_N);
        fprintf(['Leaf C balance (mmol/plant): [0, Anet, Export, Suc comsume(Rd&Nass&Leak), NSC surplus, Rd, Nass, Leak]=[%.4f, %.4f, %.4f, %.4f, '...
            '%.4f, %.4f, %.4f, %.4f]\n'],Leaf_C*12);
        fprintf(['Leaf N balance (mmol/plant): [0, N net import, N consume(Leak), NSN surplus, I-N_influx, O-N influx, O-N_outflux, Harvest(Pro), '...
            'Harvest(iN), Harvest(oN)]=[%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n'],Leaf_N);
        fprintf(['Ear C balance (mmol/plant): [0, Suc import, Suc comsume, Suc surplus, Rd, Star, Growth]=[%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, '...
            '%.4f]\n'],Ear_C*12);
        fprintf(['Ear N balance (mmol/plant): [0, N import, N consume(Growth), NSN surplus, Harvest(Pro), Harvest(oN)]=[%.4f, %.4f, %.4f, %.4f, '...
            '%.4f, %.4f]\n'],Ear_N);
        fprintf(['Stem C balance (mmol/plant): [0, Suc import, Suc comsume(Rd), Suc surplus, Harvest, Heading]=[%.4f, %.4f, %.4f, %.4f, %.4f, '...
            '%.4f]\n'],Stem_C*12);
        fprintf(['Stem N balance (mmol/plant): [0, N import, N consume(0), NSN surplus, Harvest(I-N&O-N), Harvest(Pro)]=[%.4f, %.4f, %.4f, %.4f, '...
            '%.4f, %.4f]\n'],Stem_N);
        disp('%%%%%%%%%%%%%%%%%%% end: mass balance check %%%%%%%%%%%%%%%%%%%')
        fprintf('C flow (mmol/plant): [Aleaf, Leaf_export, Root_import, Ear_import, Stem_export]=[%.4f, %.4f, %.4f, %.4f, %.4f]\n',...
            [C_abs_leaf,-Suc_in_leaf,Suc_in_root,Suc_in_ear,-Suc_in_stem]*12000);
        fprintf('Rd (mmol C/plant): [Leaf, Root, Ear, Stem]=[%.4f, %.4f, %.4f, %.4f]\n',[Rd_leaf,Rd_root,Rd_ear,Rd_stem]*12000);
        fprintf('N flow (mmol/plant): [Uroot, Leaf_export, Root_export, Ear_import, Stem_export]=[%.4f, %.4f, %.4f, %.4f, %.4f]\n',...
            [N_abs_root,-N_in_leaf,-N_shoot_root_netFlux,N_in_ear,-N_in_stem]*1000);
        fprintf('C at flowering (mmol/plant): [Leaf, Root, Ear, Stem]=[%.4f, %.4f, %.4f, %.4f]\n',...
            [NSC_leaf_flowering,NSC_root_flowering,NSC_ear_flowering,NSC_stem_flowering]*12000);
        fprintf('C at harvest (mmol/plant): [Leaf, Root, Ear, Stem]=[%.4f, %.4f, %.4f, %.4f]\n',...
            [NSC_leaf_harvest,NSC_root_harvest,NSC_ear_harvest,NSC_stem_harvest]*12000);
        fprintf('N at flowering (mmol/plant): [Leaf, Root, Ear, Stem]=[%.4f, %.4f, %.4f, %.4f]\n',...
            [N_leaf_flowering,N_root_flowering,N_ear_flowering,N_stem_flowering]*1000);
        fprintf('N at harvest (mmol/plant): [Leaf, Root, Ear, Stem]=[%.4f, %.4f, %.4f, %.4f]\n',...
            [N_leaf_harvest,N_root_harvest,N_ear_harvest,N_stem_harvest]*1000);
        fprintf(['Other N budget (mmol/plant): [Root_ass, Root_oN_import, Root_oN_export, Root_iN_export, Leaf_oN_export, Leaf_iN_import, Leaf_ass, '...
            'Leaf_oN_import, Leaf_N_sene_leak, Leaf_C_sene_leak, oN_X2P_transfer]=[%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n'],...
            other_N_budget);

        %Yield=simulation(end,516)
        Vol_solid_single_grain=(simulation(:,205)/grainNum/tillerNum/pho_Pro+simulation(:,206)/grainNum/tillerNum/pho_Star)*10^-6; % m3
        Struc_Biomass=simulation(:,544)-simulation(:,206)-simulation(:,205); % 
        Total_Biomass_solid_soluble=simulation(:,544)+(simulation(:,201)-Vol_solid_single_grain)*grainNum*tillerNum.*(simulation(:,203)*Suc_molWeight...
            + simulation(:,204)*Pro_molWeight);
        soluble_sugar=(simulation(end,201)-Vol_solid_single_grain(end))*grainNum*tillerNum*simulation(end,203)*Suc_molWeight; % g Suc
        soluble_oN=(simulation(end,201)-Vol_solid_single_grain(end))*grainNum*tillerNum*simulation(end,204)*N_molWeight; % g N
        struc_N=Struc_Biomass*grain_struc_N*N_to_Pro_mass;
        fprintf(['Grain weight at harvest (g/plant): [Total, Strucuture, Star, Pro, Sugar, AA, Structural N(Pro)]=[%.4f, %.4f, %.4f, %.4f, '...
            '%.4f, %.4f, %.4f]\n'],[Total_Biomass_solid_soluble(end),Struc_Biomass(end),simulation(end,206),simulation(end,205),soluble_sugar,...
            soluble_oN,struc_N(end)]);
        GNC=N_ear_harvest*14/Total_Biomass_solid_soluble(end)*100;
        fprintf('Grain nitrogen concentration at harvest (%%): %.4f\n',GNC);
    end
    toc
end
