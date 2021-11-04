function LRT=SCOPE_GSAT(input_var, count)

    addpath src/RTMs
    addpath src/supporting
    addpath src/fluxes
    addpath src/IO

    %% 1. define constants
    constants = define_constants;

    %% 2. paths
    path_input      = 'input/';          % path of all inputs
    path_of_code    = cd;

    %% 3. simulation options
    fid = fopen('set_parameter_filenames.csv','r');
    parameter_file = textscan(fid,'%s','Delimiter', ',');
    fclose(fid);

    fid = fopen([path_input parameter_file{1}{1}],'r');
    Ni = textscan(fid,'%d%s','Delimiter',',');%,'Whitespace','');
    fclose(fid);
    N = double(Ni{1});

    options.lite                = N(1);    % lite version
    options.calc_fluor          = N(2);    % calculate chlorophyll fluorescence in observation direction
    options.calc_planck         = N(3);    % calculate spectrum of thermal radiation
    options.calc_xanthophyllabs = N(4);    % include simulation of reflectance dependence on de-epoxydation state
    options.soilspectrum        = N(5);    % 0: use soil reflectance from file; 1: calculate soil reflectance with BSM
    options.Fluorescence_model  = N(6);     %0: empirical, with sustained NPQ (fit to Flexas' data); 1: empirical, with sigmoid for Kn; 2: Magnani 2012 model
    options.apply_T_corr        = N(7);     % correct Vcmax and rate constants for temperature in biochemical.m
    options.verify              = N(8);
    options.saveCSV             = N(9);
    options.mSCOPE              = N(10);
    options.simulation          = N(11);    % 0: individual runs (specify all input in a this file)
    % 1: time series (uses text files with meteo input as time series)
    % 2: Lookup-Table (specify the values to be included)
    % 3: Lookup-Table with random input (specify the ranges of values)
    options.calc_directional     = N(12);    % 0: calculate full BRDF (many angles)
    options.calc_vert_profiles   = N(13);
    options.soil_heat_method     = N(14);  % 0 - GAM=Soil_Inertia0(lambdas), 1 - GAM=Soil_Inertia1(SMC), 2 - G=0.35*Rn (always in no TS)
    options.calc_rss_rbs          = N(15);  % 0 - fixed, 1 calc

    if options.simulation>2 || options.simulation<0, fprintf('\n simulation option should be between 0 and 2 \r'); return, end

    %% 3. file names
    f_names = {'Simulation_Name','soil_file','optipar_file','atmos_file', 'Dataset_dir',...
        'meteo_ec_csv', 'vegetation_retrieved_csv', 'LIDF_file', 'verification_dir', ...
        'mSCOPE_csv', 'nly'};  % must be in this order
    cols = {'t', 'year', 'Rin','Rli', 'p','Ta','ea','u','RH', 'VPD', 'tts','tto', 'psi' ...  % expected from EC file as well as ('Ca','SMC')
        'Cab','Cca','Cdm','Cw','Cs','Cant','N',...  % leaf
        'SMC','BSMBrightness', 'BSMlat', 'BSMlon',...  % soil
        'LAI', 'hc', 'LIDFa', 'LIDFb',...  % canopy
        'z','Ca', ...  % meteo
        'Vcmo', 'm',...  % biochemistry;
        'atmos_names'
        };

    fnc = [f_names, cols];
    F = struct('FileID', fnc);

    fid = fopen([path_input parameter_file{1}{2}],'r');
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            charline = char(line);
            if ~(charline(1) == '%')
                X = textscan(line,'%s%s','Delimiter', ',', 'Whitespace','\t');
                x = X{1}; y = X{2};
                k = find(strcmp(fnc,x{:}));
                if ~isempty(k) && ~isempty(y)
                    F(k).FileName = y{:};
                end
            end
        end
    end
    fclose(fid);

    %% 4. input data
    k = 1;
    fid = fopen([path_input parameter_file{1}{3}], 'r');
    clear('X')
    while ~feof(fid)
        line    = fgetl(fid);
        y       = textscan(line,'%s', 'Delimiter', ',', 'TreatAsEmpty',  ' ');
        varnames(k)  = y{1}(1); %#ok<SAGROW>
        X(k).Val   = str2double(y{:});
        k       = k+1;
    end
    fclose(fid);
    V                           = assignvarnames();

    for i = 1:length(V)
        j = find(strcmp(varnames,V(i).Name));
        if isempty(j)
            if i==2
                fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data...');
                fprintf(1,'%s %s %s\n', 'I will use 0.25*Cab instead');
                options.Cca_function_of_Cab = 1;
            else
                if ~(options.simulation==1) && (i==30 || i==32)
                    fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data...');
                    fprintf(1,'%s %s %s\n', 'I will use the MODTRAN spectrum as it is');
                else
                    if (options.simulation == 1 || (~options.simulation && (i<46 || i>50 ))) && i<65
                        fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data');
                        if (options.simulation ==1 && (i==1 ||i==9||i==22||i==23||i==54 || (i>29 && i<37)))
                            fprintf(1,'%s %s %s\n', 'I will look for the values in Dataset Directory "',F(5).FileName,'"');
                        else
                            if (i== 24 || i==25)
                                fprintf(1,'%s %s %s\n', 'will estimate it from LAI, CR, CD1, Psicor, and CSSOIL');
                                options.calc_zo = 1;
                            else
                                if (i>38 && i<44)
                                    fprintf(1,'%s %s %s\n', 'will use the provided zo and d');
                                    options.calc_zo = 0;
                                else
                                    if ~((options.simulation ==1 && (i==30 ||i==32)))
                                        fprintf(1,'%s \n', 'this input is required: SCOPE ends');
                                        return
                                    elseif (options.simulation ==1 && (i==30 ||i==32))
                                        fprintf(1,'%s %s %s\n', '... no problem, I will find it in Dataset Directory "',F(5).FileName, '"');
                                    end
                                end
                            end
                        end
                    end
                end
            end
        else
            k = find(~isnan(X(j).Val));
            if ~isempty(k)
                V(i).Val = X(j).Val(k);
            else
                V(i).Val            = -999;
            end
        end
    end
    
    %% 5. revise input data
    options.verify              = 0;
    options.simulation          = 0;

    % soil
    V(17).Val     = input_var(1);      % rss->Soil resistance for evaporation
    V(54).Val   = input_var(2);      % SMC->Volumetric soil moisture content
    % leaf optical
    V(6).Val   = input_var(3);      % N->Mesophyll structural parameter in prospect
    V(4).Val  = input_var(4);       % Cw->Water content in PROSPECT
    V(3).Val= input_var(5);       % Cdm->Dry matter content in PROSPECT
    V(5).Val    = input_var(6);      % Cs->Senescence factor in PROSPECT
    V(1).Val    = input_var(7);      % Cab->Chlorophyll content in PROSPECT
    % leaf physiology(TB12-D, TB12, MD12)
    V(10).Val     = input_var(8);      % m->Ball?Berry stomatal conductance parameter
    V(12).Val   = input_var(9);      % kV->Extinction coefficient for a vertical profile of Vcmo (maximum value of Vcmo occurs at the top of the canopy).
    V(13).Val   = input_var(10);     % Rdparam->Parameter for dark respiration (Rd = Rdparam ? Vcmo)
    V(9).Val  = input_var(11);      % Vcmo->maximum carboxylation capacity (at optimum temperature)
    % Leaf physiology (MD12 only)
    V(57).Val= input_var(12);      % kNPQs->Rate constant of sustained thermal dissipation
    V(58).Val    = input_var(13);     % qLs->Fraction of functional reaction centers
    % Aerodynamic
    V(44).Val    = input_var(14);     % rbs->Soil boundary layer resistance
    V(45).Val     = input_var(15);     % rwc->Within canopy layer resistance
    V(38).Val   = input_var(16);     % rb->Leaf boundary resistance
    % Micrometeorologic
    V(33).Val   = input_var(17);     % p->Air pressure
    V(35).Val  = input_var(18);      % u->Wind speed
    V(37).Val= input_var(19);      % Oa->O2 concentration in the air
    V(34).Val    = input_var(20);     % ea->Atmospheric vapor pressure
    V(36).Val    = input_var(21);     % Ca->CO2 concentration in the air
    V(31).Val     = input_var(22);     % Ta->Air temperature
    V(30).Val   = input_var(23);     % Rin->Incoming shortwave radiation
    V(32).Val   = input_var(24);     % Rli->Incoming longwave radiation
    % Canopy
    V(28).Val   = input_var(25);     % lw->Leaf width
    V(26).Val  = input_var(26);      % LIDFa->LIDF parameter a, which controls the average leaf slope 
    V(27).Val= input_var(27);      % LIDFb->LIDF parameter b, which controls the distribution's bimodality
    V(22).Val    = input_var(28);     % LAI->Leaf area index
    V(23).Val    = input_var(29);     % hc->Canopy height
    % Geometry
    V(52).Val     = input_var(30);     % VZA->Viewing zenith angle
    V(51).Val   = input_var(31);     % SZA->Sun zenith angle
    V(53).Val   = input_var(32);     % RAA->Relative azimuth angle

    
    

    %% 6. Load spectral data for leaf and soil
    load([path_input,'fluspect_parameters/', F(3).FileName]);
    if options.soilspectrum ==0
        rsfile = load([path_input,'soil_spectra/', F(2).FileName]);        % file with soil reflectance spectra
    end

    %% 8. Define canopy structure and other 'fixed' parameters
    canopy.nlincl   = 13;
    canopy.nlazi    = 36;
    canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
    canopy.lazitab  = ( 5:10:355 );           % a row
    soilemp.SMC     = 25;        % empirical parameter (fixed) for BSM
    soilemp.film    = 0.015;     % empirical parameter (fixed) for BMS
    LIDF_file = F(8).FileName;
    if  ~isempty(LIDF_file)
        canopy.lidf     = dlmread([path_input,'leafangles/',LIDF_file],'',3,0);
    end

    %% 10. Define spectral regions
    [spectral] = define_bands;

    %% 11. load time series data
    if options.simulation == 1
        vi = ones(length(V),1);
        [soil,leafbio,canopy,meteo,angles,xyt]  = select_input(V,vi,canopy,options,constants);
        [V, xyt, mly_ts, atmo_paths]  = load_timeseries(V, F, xyt, path_input);
    else
        soil = struct;
    end

    %% 12. preparations
    %% soil heat
    if options.simulation==1
        if options.soil_heat_method<2
            if (isempty(meteo.Ta) || meteo.Ta<-273), meteo.Ta = 20; end
            soil.Tsold = meteo.Ta*ones(12,2);
        end
    end

    %% variables
    nvars = length(V);
    vmax = cellfun(@length, {V.Val})';
    vmax([14,27],1) = 1; % these are Tparam and LIDFb
    vi      = ones(nvars,1);
    switch options.simulation
        case 0, telmax = max(vmax);  [xyt.t,xyt.year]= deal(zeros(telmax,1));
        case 1, telmax = size(xyt.t,1);
        case 2, telmax  = prod(double(vmax)); [xyt.t,xyt.year]= deal(zeros(telmax,1));
    end
    % [rad,thermal,fluxes] = initialize_output_structures(spectral);

    if options.calc_directional
        anglesfile          = load([path_input,'directional/brdf_angles2.dat']); %     Multiple observation angles in case of BRDF calculation
        directional.tto     = anglesfile(:,1);              % [deg]             Observation zenith Angles for calcbrdf
        directional.psi     = anglesfile(:,2);              % [deg]             Observation zenith Angles for calcbrdf
        directional.noa     = length(directional.tto);      %                   Number of Observation Angles
    else
        directional = NaN;
    end

    %% irradiance
    atmfile = fullfile(path_input, 'radiationdata', F(4).FileName);
    if options.simulation == 1 && ~isempty(atmo_paths)
        atmfile = atmo_paths{1};
    end
    atmo = load_atmo(atmfile, spectral.SCOPEspec);

    %% 13. create output files
    [Output_dir, f, fnames] = create_output_files_binary(parameter_file, F, path_of_code, path_input, spectral,options);

    %% 14. Run the models
    fprintf('\n The calculations start now \r')
    calculate = 1;
    tic

    for k = 1:telmax
        if options.simulation == 1, vi(vmax>1) = k; end
        if options.simulation == 0, vi(vmax==telmax) = k; end
        [soil,leafbio,canopy,meteo,angles,xyt] = select_input(V,vi,canopy,options,constants,xyt,soil);
        canopy.nlayers  = ceil(10*canopy.LAI/canopy.Cv);
        canopy.nlayers = max(2, canopy.nlayers);  % patch for LAI < 0.1
        nl              = canopy.nlayers;
        x        = (-1/nl : -1/nl : -1)';         % a column vector
        canopy.xl       = [0; x];                 % add top level
        % canopy.xl(1:end-1) = canopy.xl(1:end-1)+canopy.xl(1:end-1)-1/(2*nl); % middle of the thin layer

        if options.simulation ~=1
            fprintf('simulation %i ', k );
            fprintf('of %i \n', telmax);
        else
            calculate = ~isnan(meteo.p*meteo.Ta*meteo.ea*meteo.u.*meteo.Rin.*meteo.Rli);
            fprintf('time = %s: %i / %i\n', datestr(xyt.t(k)), k, telmax)
            if isnan(meteo.p*meteo.Ta*meteo.ea*meteo.u.*meteo.Rin.*meteo.Rli)
                warning('run is invalid: there is NaN somewhere in meteo input [p, Ta, ea, u, Rin, Rli]')
            end
        end

        if calculate
            if isempty(LIDF_file)
                canopy.lidf     = leafangles(canopy.LIDFa,canopy.LIDFb);    % This is 'ladgen' in the original SAIL model,
            end

            %% leaf radiative transfer model FLUSPECT
            leafbio.emis        = 1-leafbio.rho_thermal-leafbio.tau_thermal;
            leafbio.V2Z         = 0;

            if options.simulation == 1 && ~isempty(fieldnames(mly_ts))  % means that options.simulation == 1
                mly.nly    = mly_ts.nly;
                mly.pLAI   = mly_ts.pLAI(k, :);
                mly.totLAI = sum(mly.pLAI);
                mly.pCab   = mly_ts.pCab(k, :);
                mly.pCca   = mly_ts.pCca(k, :);
                mly.pCdm   = mly_ts.pCw(k, :);
                mly.pCw    = mly_ts.pCw(k, :);
                mly.pCs    = mly_ts.pCs(k, :);
                mly.pN     = mly_ts.pN(k, :);
            elseif k == 1 && options.mSCOPE
                mly = input_mSCOPE(fullfile('input', 'mSCOPE.csv'));
            else
                if options.mSCOPE
                    warning('I do not know how to use mSCOPE layers with multiple but non time series runs, so I will not use it')
                end
                mly.nly      = 1;
                mly.pLAI     = canopy.LAI;
                mly.totLAI   = canopy.LAI;
                mly.pCab     = leafbio.Cab;
                mly.pCca     = leafbio.Cca;
                mly.pCdm     = leafbio.Cdm;
                mly.pCw      = leafbio.Cw;
                mly.pCs      = leafbio.Cs;
                mly.pN       = leafbio.N;
            end

            if options.simulation == 1 && ~isempty(atmo_paths) && k > 1
                atmfile_k = atmo_paths{k};
                if ~strcmp(atmfile_k, atmo_paths{k-1})
                    atmo = load_atmo(atmfile_k, spectral.SCOPEspec);
                end
            end

            leafopt = fluspect_mSCOPE(mly,spectral,leafbio,optipar, nl);
            leafopt.refl(:, spectral.IwlT) = leafbio.rho_thermal;
            leafopt.tran(:, spectral.IwlT) = leafbio.tau_thermal;

            if options.calc_xanthophyllabs
                leafbio.V2Z     = 1;
                leafoptZ = fluspect_mSCOPE(mly,spectral,leafbio,optipar, nl);
                leafopt.reflZ   = leafopt.refl;
                leafopt.tranZ   = leafopt.tran;
                leafopt.reflZ(:, spectral.IwlP) = leafoptZ.refl(:, spectral.IwlP);
                leafopt.tranZ(:, spectral.IwlP) = leafoptZ.tran(:, spectral.IwlP);
            end

            %% soil reflectance model BSM
            if options.soilspectrum == 0
                soil.refl       = rsfile(:,soil.spectrum+1);
            else
                soil.refl       = BSM(soil,optipar,soilemp);
            end
            soil.refl(spectral.IwlT) = soil.rs_thermal;

            %% four stream canopy radiative transfer model for incident radiation
            [rad,gap,profiles]       = RTMo(spectral,atmo,soil,leafopt,canopy,angles,constants,meteo,options);

            %% energy balance
            [iter,rad,thermal,soil,bcu,bch,fluxes]             ...
                = ebal(constants,options,rad,gap,  ...
                meteo,soil,canopy,leafbio, k, xyt);

            %         [iter,rad,thermal,soil,bcu,bch,fluxes]             ...
            %             = ebal_sunshade(constants,options,rad,gap,  ...
            %             meteo,soil,canopy,leafbio);

            %         [iter,rad,thermal,soil,bcu,bch,fluxes]             ...
            %             = ebal_bigleaf(constants,options,rad,gap,  ...
            %             meteo,soil,canopy,leafbio);

            %% fluorescence radiative transfer model
            if options.calc_fluor
                [rad]           = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles,bcu.eta,bch.eta);
            end

            %% radiative transfer model for PRI effects
            if options.calc_xanthophyllabs
                [rad] = RTMz(constants,spectral,rad,soil,leafopt,canopy,gap,angles,bcu.Kn,bch.Kn);
            end

            rad  = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,thermal.Tcu,thermal.Tch,thermal.Tsu,thermal.Tsh,1,spectral);
            if options.calc_planck
                rad = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,thermal.Tcu,thermal.Tch,thermal.Tsu,thermal.Tsh);
            end

            %% computation of data products
            % aPAR, LST, NPQ, ETR, photosynthesis, SIF-reabsorption correction

            % aPAR [umol m-2 s-1, total canopy and total chlorphyll]
            switch options.lite
                case 0, integr = 'angles_and_layers';
                otherwise, integr = 'layers';
            end

            Ps = gap.Ps(1:nl);
            Ph = (1-Ps);

            canopy.Pntot_Cab = canopy.LAI*(meanleaf(canopy,rad.Pnh_Cab,'layers',Ph)+meanleaf(canopy,rad.Pnu_Cab,integr,Ps)); % net PAR Cab leaves
            canopy.Pntot     = canopy.LAI*(meanleaf(canopy,rad.Pnh,'layers',Ph)+meanleaf(canopy,rad.Pnu,integr,Ps)); % net PAR leaves
            canopy.Rntot_Cab = canopy.LAI*(meanleaf(canopy,rad.Rnh_Cab,'layers',Ph)+meanleaf(canopy,rad.Rnu_Cab,integr,Ps)); % net PAR Cab leaves
            canopy.Rntot_PAR = canopy.LAI*(meanleaf(canopy,rad.Rnh_PAR,'layers',Ph)+meanleaf(canopy,rad.Rnu_PAR,integr,Ps)); % net PAR Cab leaves

            % LST [K] (directional, but assuming black-body surface!)
            canopy.LST      = (pi*(rad.Lot+rad.Lote)./(constants.sigmaSB)).^0.25;

            % photosynthesis [mumol m-2 s-1]
            canopy.A        = canopy.LAI*(meanleaf(canopy,bch.A,'layers',Ph)+meanleaf(canopy,bcu.A,integr,Ps)); % photosynthesis

            % electron transport rate [mumol m-2 s-1]
            canopy.Ja       = canopy.LAI*(meanleaf(canopy,bch.Ja,'layers',Ph)+meanleaf(canopy,bcu.Ja,integr,Ps)); % electron transport

            % non-photochemical quenching (energy) [W m-2]
            canopy.ENPQ     = canopy.LAI*(meanleaf(canopy,rad.Rnh_Cab.*bch.Phi_N,'layers',Ph)+meanleaf(canopy,rad.Rnu_Cab.*bcu.Phi_N,integr,Ps)); % NPQ energy;

            % computation of re-absorption corrected fluorescence
            % Yang and Van der Tol (2019); Van der Tol et al. (2019)
            aPAR_Cab_eta    = canopy.LAI*(meanleaf(canopy,bch.eta .* rad.Rnh_Cab,'layers',Ph)+meanleaf(canopy,bcu.eta .* rad.Rnu_Cab,integr,Ps)); %
            if options.calc_fluor
                rad.EoutFrc     = leafbio.fqe*aPAR_Cab_eta;
                rad.EoutFrc_    = 1E3*rad.EoutFrc*optipar.phi(spectral.IwlF);
                rad.sigmaF      = pi*rad.LoF_./rad.EoutFrc_;
            end

            rad.Lotot_      = rad.Lo_+rad.Lot_;
            rad.Eout_       = rad.Eout_+rad.Eoutte_;
            if options.calc_fluor
                rad.Lototf_     = rad.Lotot_;
                rad.Lototf_(spectral.IwlF') = rad.Lototf_(spectral.IwlF)+rad.LoF_;
            end

            if options.calc_directional
                directional = calc_brdf(constants,options,directional,spectral,angles,atmo,soil,leafopt,canopy,meteo,thermal,bcu,bch);
                savebrdfoutput(options,directional,angles,spectral,Output_dir)
            end
            % reflectance
            canopy.reflectance     = pi*rad.Lo_./(rad.Esun_+rad.Esky_);

            rad.Lo = 0.001 * Sint(rad.Lo_(spectral.IwlP),spectral.wlP);
            %% write output
            n_col = output_data_binary(f, k, xyt, rad, canopy, V, vi, vmax, options, fluxes, meteo, iter);

            %% update input
            if options.simulation==2 && telmax>1, vi  = count_k(nvars,vi,vmax,1); end

        end
    end
    toc

LRT     = rad.LoF_(count);
