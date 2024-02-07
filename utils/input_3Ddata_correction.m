function [saturation_obs, porosity_prof] = input_3Ddata_correction(model, saturation_obs, porosity_prof)
% INPUT_3DDATA_CORRECTION Corrects 3D saturation and porosity data based on grid dimensions.
%
% DESCRIPTION:
%   This function adjusts the input 3D saturation observations and porosity profiles to align with the modified
%   grid dimensions of a model. It ensures that the saturation and porosity data fit within the new grid
%   dimensions by resizing the data arrays if necessary and removing data corresponding to boundary cells.
%
% SYNOPSIS:
%   [saturation_obs, porosity_prof] = input_3Ddata_correction(model, saturation_obs, porosity_prof)
%
% PARAMETERS:
%   model - Struct containing simulation and grid information, modified to include new grid dimensions.
%   saturation_obs - Struct containing the saturation observations data:
%                    - time: Array of observation times
%                    - saturation_array: Cell array of saturation profiles for each time step
%                    - saturation_matrix: Cell array of 3D saturation matrices for each time step
%   porosity_prof - Struct containing porosity profiles data:
%                   - porosity_profile_3d: 3D array of porosity data
%                   - porosity_profile_vertical: Vector of porosity data reshaped from 3D array
%
% RETURNS:
%   saturation_obs - Updated struct with corrected saturation data to fit the model's grid dimensions.
%   porosity_prof - Updated struct with corrected porosity profiles to fit the model's grid dimensions.
%
% EXAMPLES:
%   Assuming `model`, `saturation_obs`, and `porosity_prof` are already defined and `model` has been updated
%   to include a new grid via `CreateGrid` function:
%   [saturation_obs, porosity_prof] = input_3Ddata_correction(model, saturation_obs, porosity_prof);
%
% NOTES:
%   - This function assumes that the input model struct has been previously updated to include a modified
%     grid structure, typically by using the CreateGrid function.
%   - The correction process involves resizing the 3D matrices of saturation and porosity data to match the
%     new grid dimensions and removing entries that correspond to boundary cells that are not part of the
%     simulation grid anymore.
%
%
% ----------------------------------
% (c) 2024
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%

% pre-check the saturation observation to be between 0 and 1
% not a good idea since it messes your average slice saturations
% eps = 0;
% for i = 1:length(saturation_obs.time)
%     saturation_obs.saturation_array{1, i}(saturation_obs.saturation_array{1, i} < 0) = 0 + eps;
%     saturation_obs.saturation_array{1, i}(saturation_obs.saturation_array{1, i} > 1) = 1 - eps;
%     saturation_obs.saturation_matrix{1, i}(saturation_obs.saturation_matrix{1, i} < 0) = 0 + eps;
%     saturation_obs.saturation_matrix{1, i}(saturation_obs.saturation_matrix{1, i} > 1) = 1 - eps;
% end

t_obs = saturation_obs.time;
z_matrix_vertical = saturation_obs.saturation_array;
z_matrix_3d_saved = saturation_obs.saturation_matrix;
porosity_profile_3d = porosity_prof.porosity_profile_3d;
porosity_profile_vertical = porosity_prof.porosity_profile_vertical;

% 2 is coming from number of boundary cells
model_copy = model;
model_copy.simulation.nCells.value = model_copy.simulation.nCells.value - 2;
model_copy = Create3DGrid(model_copy);
remove_mask = model_copy.grid.remove_mask;

grid_dims = model.grid.G.cartDims;
% because of boundary cells
grid_dims(1) = grid_dims(1) - 2;
if grid_dims(1) < length(z_matrix_3d_saved{1, 1}) || grid_dims(2) < width(z_matrix_3d_saved{1, 1})
    for i = 1:length(t_obs)
        z_matrix_3d_saved{1, i} = imresize3(z_matrix_3d_saved{1, i}, grid_dims);
        z_matrix_vertical{1, i} = z_matrix_3d_saved{1, i}(:);
        z_matrix_vertical{1, i}(remove_mask) = [];
    end
else
    for i = 1:length(t_obs)
        z_matrix_vertical{1, i}(remove_mask) = [];
    end
end

if grid_dims(1) < length(porosity_profile_3d) || grid_dims(2) < width(porosity_profile_3d)
    porosity_profile_3d = imresize3(porosity_profile_3d, grid_dims);
    porosity_profile_vertical = porosity_profile_3d(:);
    porosity_profile_vertical(remove_mask) = [];
else
    porosity_profile_vertical(remove_mask) = [];
end




saturation_obs.time = t_obs;
saturation_obs.saturation_array = z_matrix_vertical;
saturation_obs.saturation_matrix = z_matrix_3d_saved;
porosity_prof.porosity_profile_3d = porosity_profile_3d;
porosity_prof.porosity_profile_vertical = porosity_profile_vertical;