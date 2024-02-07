%% Simulation of 3d heterogenous cores
%
% DESCRIPTION: Example of how to run a case with the SCAL module
% Have a look at "docs/list of examples.docx" for a list of available input
% examples + "doc/List of keywords in the settings files.docx" for a
% documentation about what is the role of each keyword
% ----------------------------------
% (c) 2024
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%% clear memory, close all figures and screen
clc;
clear all;
close all;
clear classes;

%% configure model from file
% adding mrst modules
mrstModule add ad-core ad-props ad-blackoil ad-scal-heterogeneous mrst-gui
% input absolute path to the settings file containing the keywords and
% input parameters/tables
model = Configure(fullfile (ROOTDIR,'modules','ad-scal-heterogeneous','examples','SettingsDecaneBrine3D.txt'));

%% App function
model.App.include = false;

%% verbose to display internal messages
% set true to show detailed messages about the simulation
verbose = true; 
model.verbose = verbose;

%% discretize geometry
% create the grid on which the simulation is run
model = Create3DGrid(model);

%% create rock model
% setting porosity and permeability data for the simulation
model = Create3DRock(model);

%% create saturation functions model
% creating the capillary pressure and relative permeability tables for the
% simulation
model = CreatePc(model);
model = CreateKr(model);

%% plot satuartion functions and comparison plots
% plot the pc and kr tables created in the previous step for verification
% and comparison purposes 
if not(strcmpi(model.simulation.type,strcat('historymatch')))
    PlotKrPc(model);
    if or(isfield(model.experiment.satfun,'kr_compare'), ...
            isfield(model.experiment.satfun,'pc_compare'))
        PlotKrPcComparison(model)
    end
    % just a pause before you run the simulations to check if you have the
    % correct input relative permeability and capillary pressures for 
    % your simulations - you can also plot them against other kr and pc 
    % curves for comparison purposes
    % disp('Check input saturation functions, continue? [enter for yes/control + c to no]');
    % pause()
end

%% create fluid model
% create functions to interpolate the pc and kr saturation function tables
model = CreateFluid(model);

%% run Froward/Histormatch simulation
switch lower(model.simulation.type)
    case 'forward'   
% ------------forward heterogenous model-----------------------------------
        if model.experiment.rock.heterogeneous
            fprintf('Simulating the prior homogeneous model...\n')
            model = Run3D(model, false);
            fprintf('calculating the pc scaling factors...\n')
            alpha = model.experiment.rock.alpha.value;
            matFolder_Decane = '.\decane_brine\';
            matFolder_CO2 = '.\co2_brine\';
            load(fullfile (ROOTDIR, 'modules','ad-scal-heterogeneous','examples', matFolder_Decane, 'time_array.mat'))
            t_obs = time_array;
            load(fullfile (ROOTDIR, 'modules','ad-scal-heterogeneous','examples',matFolder_Decane, 'z_matrix_3d_saved_satprof.mat'))
            load(fullfile (ROOTDIR, 'modules','ad-scal-heterogeneous','examples',matFolder_Decane, 'z_matrix_vertical_satprof.mat'))
            saturation_obs.time = t_obs;
            saturation_obs.saturation_array = z_matrix_vertical;
            saturation_obs.saturation_matrix = z_matrix_3d_saved;
            load(fullfile (ROOTDIR, 'modules','ad-scal-heterogeneous','examples', matFolder_Decane, 'porosity_profile_3d.mat'))
            load(fullfile (ROOTDIR, 'modules','ad-scal-heterogeneous','examples', matFolder_Decane, 'porosity_profile_vertical.mat'))
            porosity_prof.porosity_profile_3d = porosity_profile_3d;
            porosity_prof.porosity_profile_vertical = porosity_profile_vertical;
            [saturation_obs, porosity_prof] = input_3Ddata_correction(model, saturation_obs, porosity_prof);
            het_index_mask = model.experiment.rock.het_index_mask;
            f_factor = f_factor_calculator_3d(model, alpha, true, ...
                het_index_mask, saturation_obs, porosity_prof);
            % sw_diff = calculate_sw_dif(model, f_factor, saturation_obs, het_index_mask);
            % sw_diff = zeros(size(f_factor));
            % disp("Check f factor calculations, continue? [enter for yes/control + c to no]");
            % pause()

            model = Create3DGridHeterogeneous(model);
            model = Create3DRockHeterogeneous(model, porosity_prof, f_factor);
            model = Create3DFluidHeterogenous(model, f_factor, sw_diff);
            fprintf('Simulating the heterogeneous model...\n')
            model = Run3D(model, true);
%             disp("Plot saturation profile? [enter for yes/ control + c for no]");
%             pause()
            Plot3DSaturationProfiles(model, het_index_mask, saturation_obs)
            Compare3DSaturationProfiles(model, het_index_mask, saturation_obs)
        else
% ------------forward homogeneous model-----------------------------------
            fprintf('Simulating the homogeneous model...\n')
            model = Run3D(model, true);
%             choice_4 = input("Plot saturation profile? [y/n]", "s");
%             if strcmpi(choice_4, "y")
                plot_saturation_profile(model)
%             end
        end
    case 'historymatch'
% ------------investigate the input range for------------------------------
%     plot_input_range(model)
%     fprintf('Want to continue? Press Enter to continue \n')
%     pause
%-------------Exploratory data analysis------------------------------------
        filter_no.pressure = 50;
        filter_no.prod = 10;
        filtered_data = PlotObservation_pre_hm(model, filter_no);
        plot_saturation_profile_pre_hm(model)
        choice = input('History match with the filtered data? [y/n]','s');
        if strcmpi(choice, 'y')
            model = replace_observation_with_filtered_data(model, filtered_data);
        end
        close
%--------------------------------------------------------------------------
        fprintf('Starting the history match...\n')
%----------Multi objective history matching ------------------------------
        switch model.history_match.algorithm
            case 'ga_multi_objective'
                [model, model_cent] = historymatch_multi_obj(model);
                [x_best, fval_sum_best] = get_x_best(model);
                model.history_match.x = x_best;
                model.history_match.fval = fval_sum_best;
            otherwise
%----------single objective history matching ------------------------------
        [model, model_cent] = historymatch(model);
        end
%% plotting after history match
    plot_HM_results(model)
    obj_fun = model.history_match.obj_fun;
    switch lower(obj_fun)
        case 'simultaneous'
            plot_HM_response(model,'ss')
            model_cent.history_match = model.history_match;
            plot_HM_response(model_cent,'cent')
        otherwise
            if model.experiment.rock.heterogeneous
                plot_HM_response(model,'', model.experiment.rock.het_index_mask)
            else
                plot_HM_response(model,'')
            end
    end
end

%% output results
% if not(isempty(model.output))
%     if(model.output.include)
%         SaveResults(model);
%     end
% end

%% compute and plot fractional flow
% [sw, fw] = compute_fractional_flow (model);
% plot_fractional_flow(sw, fw);
