function Plot3DSaturationProfiles(model, index_mask, saturation_obs)
%
% DESCRIPTION: 
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
t_obs = saturation_obs.time*60;

% get the saturation from the simulation
t = model.history_match.Sw_profile(2:end,1);

mutual_idx_t = ismember(t,t_obs, 'row');
mutual_idx_t_obs = ismember(t_obs,t, 'row');

saturation_3D_matrix = saturation_obs.saturation_matrix;
equivalent_3D_obs_matrix = saturation_3D_matrix(mutual_idx_t_obs);
indexed_3D_obs = equivalent_3D_obs_matrix(index_mask);
new_mutual_idx_t = [0 ;mutual_idx_t];
indexed_slide_average_model = model.history_match.Sw_profile(logical(new_mutual_idx_t),:);
indexed_slide_average_model = indexed_slide_average_model(index_mask, :);

index_mask_time = t_obs(mutual_idx_t_obs) / 3600;
index_mask_time = index_mask_time(index_mask);

t_1D = model.experiment.observation.satProfile.table(2:end,1);
mutual_idx_t_obs_1D = ismember(t_1D,t_obs, 'row') & ismember(t_1D,t, 'row');
satObs_1D_mutual = model.experiment.observation.satProfile.table(logical([0; mutual_idx_t_obs_1D]), :);
satObs_1D_mutual_indexed = satObs_1D_mutual(index_mask,:);

f = figure;
f.Position = [100 100 1000 400];
G = model.grid.G;
n_slice_cells = length(G.cells.centroids(G.cells.centroids(:,1)==G.cells.centroids(3,1)));
for i = 1:length(index_mask)
    n_cells = model.grid.G.cartDims(1) - 2;
%     s(i) = plot(1:n_cells,  sum(sum(indexed_3D_obs{i},3),2)/n_slice_cells, ...
%         'DisplayName', strcat("Slice average of 3D observations at ", ...
%         string(index_mask_time(i)), " hours"));
    hold on
    p(i) = plot(1:n_cells, 1- indexed_slide_average_model(i,3:end-1), '-', 'markersize', 2, ...
        'DisplayName', 'Slice average of heterogeneous 3D simulations');
%     p(i).Color = s(i).Color;
    sc(i) = plot(linspace(1,n_cells,length(model.experiment.observation.satProfile.table)-1), 1- satObs_1D_mutual_indexed(i,2:end), 'o',  ...
        'DisplayName', '1D observations');
    p(i).Color = sc(i).Color;
end
legend('Location','bestoutside')
