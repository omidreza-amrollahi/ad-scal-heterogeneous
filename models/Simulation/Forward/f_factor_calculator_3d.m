function f_factor = f_factor_calculator_3d(model, alpha, plots, index_mask,...
    saturation_obs, porosity_profile)
%
% DESCRIPTION: Calculates the scaling factors for the capillary pressure
% heterogenity scaling
%
% SYNOPSIS:
%   f_factor = f_factor_calculator(model, alpha, plots, index_mask)
%
% PARAMETERS:
%   - model - main modeling struct of the simulation
%   - alpha - between 0 and 1 defining the objective function for the
%   scaling factor calculations, with 1 moving towards pressure and 0
%   moving towards saturation
%   - plots - activates diagnostics plots
%   - index_mask - index of the saturation profiles to use for calculations
%
% RETURNS:
%   f_factor - calculated capillary scaling factors
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%

    % get the saturation from the observation
    sw_obs = saturation_obs.saturation_array;
    t_obs = saturation_obs.time*60;
    
    % get the saturation from the simulation
    t = model.history_match.Sw_profile(2:end,1);
    states{1} = model.state0;
    states = [states; model.dynamic.states];

    mutual_idx_t = ismember(t,t_obs, 'row');
    equivalet_states = states(mutual_idx_t);
    mutual_idx_t_obs = ismember(t_obs,t, 'row');
    equivalet_so_obs = sw_obs(mutual_idx_t_obs);
    equivalet_sw_obs = cell(length(equivalet_so_obs),1);
    for i = 1:length(equivalet_so_obs)
        equivalet_sw_obs{i} = 1 - equivalet_so_obs{i};
    end
    

    sw_model_matrix = []; sw_obs_matrix = [];
    for i = index_mask
        sw_model_matrix = [sw_model_matrix, equivalet_states{i,1}.s(:,1)];
        sw_obs_matrix = [sw_obs_matrix, equivalet_sw_obs{i}];
    end
    % obs data are in decane saturation so we need the -1 here
    % sw_obs_matrix = 1 - sw_obs_matrix;
    % assume the ones that are 1 are out of the core circle so we remove
    % them, this is because of the preprocessing step of the data
%     sw_obs_matrix(sw_obs_matrix(:,1)==1,:) = [];

    inner_mask = not(model.grid.G.inlet_mask) & not(model.grid.G.outlet_mask);
    sw_model_matrix = sw_model_matrix(inner_mask, :);

    pc_model = model.fluid.pcOW{1,2}(sw_model_matrix);
    pc_sw_obs = model.fluid.pcOW{1,2}(sw_obs_matrix);

    % data correction only for f calculations - positive
    PcAt1 = model.fluid.pcOW{1,2}(1);
    pc_sw_obs(sw_obs_matrix > 1) = PcAt1;

    constant_porosity = model.experiment.rock.poro.value;
    constant_perm = model.experiment.rock.perm.value;
    poro_array = porosity_profile.porosity_profile_vertical;
    % similar assumtion of the saturation for the porosity as well
%     poro_array(poro_array==0) = [];
%% plotting using MRST tools
    if plots
        [max_x, max_y, max_z] = size(saturation_obs.saturation_matrix{1,1});
        core_length = model.experiment.geometry.length.value / Convert('cm');
        core_diameter = model.experiment.geometry.diameter.value / Convert('cm');
        a = linspace(0, core_length, max_x+1);
        b = linspace(0, core_diameter, max_y+1);
        c = linspace(0, core_diameter, max_z+1);
        G = tensorGrid(a, b, c);
        G = computeGeometry(G);
        radius = core_diameter/2;
        euclidean_distances = sqrt((G.cells.centroids(:,2)-radius).^2+(G.cells.centroids(:,3)-radius).^2);
        remove_mask = euclidean_distances > radius;
        G = removeCells(G, remove_mask);
        plotting_struct.sw_obs_matrix = num2cell(sw_obs_matrix,1);
        plotting_struct.sw_obs_matrix_full = equivalet_sw_obs;
        plotting_struct.sw_model_matrix = num2cell(sw_model_matrix,1);
        plotting_struct.pc_sw_obs = num2cell(pc_sw_obs,1);
        plotting_struct.pc_model = num2cell(pc_model,1);
        plotting_struct.porosity = num2cell(poro_array,1);
    end 
%% initialize f with a guess coming from the first pressure profile
    f0 = pc_sw_obs(:,1) ./ pc_model(:,1);
%     f0 = pc_sw_obs(:,2) ./ pc_model(:,2);
%     f0 = (f1+f2)./2;
    obj_fun = @(x) objective_function(x, pc_model, ...
        pc_sw_obs, sw_model_matrix,...
        sw_obs_matrix, alpha);
%     inequalities = @(iterator) mycon(iterator, poro_array(:), constant_porosity, constant_perm);
    options_fmincon = optimoptions('fmincon','Display','iter');
    lb = ones(length(f0),1)*0; ub = ones(length(f0),1)*500;
    f_factor = fmincon(obj_fun,f0(:),[],[],[],[],lb,ub,[],options_fmincon);

%     PSoptions = optimoptions(@patternsearch,'Display','iter', 'UseParallel', true);
%     f_factor = patternsearch(obj_fun,f0(:),[],[],[],[],lb,ub,PSoptions);
%% in case one wants to test other optimizers
%     nvars = length(f0);
%     lb = ones(nvars, 1)*0.5;
%     ub = ones(nvars, 1)*1.5;
%     options_swarm = optimoptions('particleswarm','Display','iter');
%     f_factor = particleswarm(obj_fun,nvars,lb,ub,options_swarm);
%%
    if plots
        figure
        pc_model_p = plot(model.satfun.sw_pc, model.satfun.pc / Convert('bar'),'r-', 'DisplayName', 'Pc model');
        title('Scaled capillary pressure');
        hold on
        scaled_pc =  f_factor(:) .* pc_model;
        scale_data_p = plot(sw_obs_matrix(:), scaled_pc(:) / Convert('bar'), 'bs', 'DisplayName', 'Scaled data');
        legend([pc_model_p(1) scale_data_p(1)],{'Pc model', 'Scaled data'})
        xlabel('Water Saturation'); ylabel('Capillary pressure (bar)')
        xlim([0 1.3]); set(gca, 'YScale', 'log');
        ylim([0.03 1]); xline(1);

        figure
        title('Capillary pressure before scaling');
        hold on
        scaled_pc =  f_factor(:) .* pc_model;
        unscaled_data_p = plot(sw_obs_matrix(:), pc_model(:) / Convert('bar'), 'ks', 'DisplayName', 'Unscaled data');
        legend([pc_model_p(1) unscaled_data_p(1)],{'Pc model', 'Unscaled data'})
        xlabel('Water Saturation'); ylabel('Capillary pressure (bar)')
        xlim([0 1.3]); set(gca, 'YScale', 'log');
        ylim([0.03 1]); xline(1);

        saturation_3D_matrix = saturation_obs.saturation_matrix;
        equivalent_3D_obs_matrix = saturation_3D_matrix(mutual_idx_t_obs);
        indexed_3D_obs = equivalent_3D_obs_matrix(index_mask);

        new_mutual_idx_t = [0 ;mutual_idx_t];
        indexed_slide_average_model = model.history_match.Sw_profile(logical(new_mutual_idx_t),:);
        indexed_slide_average_model = indexed_slide_average_model(index_mask, :);

        index_mask_time = t_obs(mutual_idx_t_obs) / 3600;
        index_mask_time = index_mask_time(index_mask);
        fmt = ['You are using the saturation profiles at ' repmat(' %.2f ',...
        1,numel(index_mask_time)) ' (hours) for f factor calculations\n'];
        fprintf(fmt, index_mask_time)

        t_1D = model.experiment.observation.satProfile.table(2:end,1);
        mutual_idx_t_obs_1D = ismember(t_1D,t_obs, 'row'); % & ismember(t_1D,t, 'row');
        satObs_1D_mutual = model.experiment.observation.satProfile.table(logical([0; mutual_idx_t_obs_1D]), :);
        satObs_1D_mutual_indexed = satObs_1D_mutual(index_mask,:);

        f = figure;
        f.Position = [100 100 1000 400];
        n_slice_cells = length(G.cells.centroids(G.cells.centroids(:,1)==G.cells.centroids(3,1)));
        for i = 1:length(index_mask)
            n_cells = model.grid.G.cartDims(1) - 2;
            s(i) = plot(1:n_cells, 1 - sum(sum(indexed_3D_obs{i},3),2)/n_slice_cells, ...
                'DisplayName', strcat("Slice average of 3D observations at ", ...
                string(index_mask_time(i)), " hours"));
            hold on
            p(i) = plot(1:n_cells, indexed_slide_average_model(i,3:end-1), '--', ...
                'DisplayName', 'Slice average of homogeneous 3D simulations');
            p(i).Color = s(i).Color;
            sc(i) = scatter(linspace(1,n_cells,length(model.experiment.observation.satProfile.table)-1), satObs_1D_mutual_indexed(i,2:end), ...
                'DisplayName', '1D observations');
        end
        legend('Location','bestoutside')

        figure
        f_factor_plot = f_factor;
        plotting_struct.f_factor = num2cell(f_factor_plot,1);
        plotting_struct.scaled_pc = num2cell(scaled_pc,1);
        plotToolbar(G, plotting_struct);
        xlabel('Core Length'); ylabel('y axis'); zlabel('z axis');
        grid on; colorbar
        colormap parula(20)

    end

end

function error = objective_function(f, pc_model, pc_sw_obs, sw_model, sw_obs, alpha)
    f = f(:);
    % alpha defines the weighting trade of between the pressure and the
    % saturation
    n = size(pc_model);
    f_copy = repmat(f,1,n(2));
    difference_p = alpha * (pc_sw_obs - pc_model .* f_copy) ./  pc_sw_obs;
    difference_sw = (1-alpha) * (sw_obs - sw_model) ./  sw_obs;
    difference_total = difference_p + difference_sw;
    difference_total(isnan(difference_total)|isinf(difference_total))=[];
    error = rms(difference_total, 'all');
end
function [c,ceq] = mycon(iterator, poro_array, constant_porosity, constant_perm)
    iterator = iterator(:);
    c(1) = abs(mean(iterator) - 1) - 0.3;
    perm_array = constant_perm .* iterator .^2 .* poro_array ./ constant_porosity;
    k_hm = harmmean(perm_array);
    c(2) = abs(k_hm - constant_perm) ./ constant_perm - 0.2;
    ceq =  [];
end