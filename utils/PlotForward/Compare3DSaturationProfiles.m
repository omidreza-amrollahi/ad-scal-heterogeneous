function Compare3DSaturationProfiles(model, index_mask, saturation_obs)

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
%sw_obs_matrix(sw_obs_matrix(:,1)==1,:) = [];

inner_mask = not(model.grid.G.inlet_mask) & not(model.grid.G.outlet_mask);
sw_model_matrix = sw_model_matrix(inner_mask, :);

index_mask_time = t_obs(mutual_idx_t_obs) / 3600;
index_mask_time = index_mask_time(index_mask);

for i=1:length(index_mask)
    figure
    X = [sw_model_matrix(:,i), sw_obs_matrix(:,i)];
    corrplot(X, VarNames=["Model Prediction" "Experimental Measurements"]);
    title('Compare saturation profiles at ' + string(index_mask_time(i)) + ' hours')
end
figure
tiledlayout('flow')
for i=1:length(index_mask)
    ax = nexttile;
    histogram(ax, sw_model_matrix(:,i), 'DisplayName', 'Model'); hold on
    histogram(ax, sw_obs_matrix(:,i), 'DisplayName', 'Experiment'); legend;
    title(ax, 'Compare saturation profiles at ' + string(index_mask_time(i)) + ' hours')
end



