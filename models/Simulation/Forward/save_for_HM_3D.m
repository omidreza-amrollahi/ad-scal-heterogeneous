function model = save_for_HM_3D (model, state0)
%
% DESCRIPTION: Save the required parameters, used during the history
%              matching
%
% SYNOPSIS:
%   model = save_for_HM (model, state0)
%
% PARAMETERS:
%   model - struct containing the following fields:
%   - dynamic: information saved during the simulation
%   - experiment: saturation functions used for forward modeling
%   - grid
%   state0 - initial state of the simulation
%
% RETURNS:
%   model - struct containing the following fields:
%   - history_match: parameters used in gradient based and MCMC
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
model.history_match.calculatedp = [];
params = model.dynamic.params;
G = model.grid.G; n_cells = G.cells.num;
model.history_match.endofsched_t = model.experiment.schedule.procedure{:,2};
sample_L = model.experiment.geometry.length.value;
states = model.dynamic.states;
centroids = G.cells.centroids;
x_dims = G.cartDims(1); 

if (isfield(model.simulation,'bCells'))
    % for now only implemented for the case with boundary cells

    sw_profile_3D = zeros(length(params.cumScheduleSteps) + 1, n_cells + 1);
    sw_profile_3D(2:end,1) = params.cumScheduleSteps;
    LocationVector = repmat([0;centroids(2:x_dims-1,1) - 2*centroids(1,1);sample_L], length(centroids)/x_dims, 1);
    sw_profile_3D(1,2:end) = LocationVector';
%     sw_profile_3D(2,2) = interp1(centroids(2:3,1),state0.s(2:3,1),0,'linear','extrap');
%     sw_profile_3D(2,end) = interp1(centroids(end-2:end-1,1) - ...
%         2 * centroids(1,1),state0.s(end-2:end-1,1),sample_L,'linear','extrap');
    sw_profile_3D(2,2:end) = state0.s(:,1)';
    for i = 1:length(states)
        sw_profile_3D(i + 2, 2:end) = states{i, 1}.s(:,1)';
%         sw_profile_3D(i + 2,2) = interp1(centroids(2:3,1),states{i, 1}.s(2:3,1)',0,'linear','extrap');
%         sw_profile_3D(i + 2,end) = interp1(centroids(end-2:end-1,1) - ...
%             2 * centroids(1,1),states{i, 1}.s(end-2:end-1,1),sample_L,'linear','extrap');
    end
    
    pressure_prof_3D = zeros(length(params.cumScheduleSteps) + 1, n_cells + 1);
    pressure_prof_3D(2:end,1) = params.cumScheduleSteps;
    pressure_prof_3D(1,2:end) = LocationVector';
    pressure_prof_3D(2,2) = interp1(centroids(2:3,1),state0.pressure(2:3,1),0,'linear','extrap');
    pressure_prof_3D(2,end) = interp1(centroids(end-2:end-1,1) - ...
        2 * centroids(1,1),state0.pressure(end-2:end-1,1),sample_L,'linear','extrap');
    pressure_prof_3D(2,3:end-1) = state0.pressure(2:end-1,1)';
    for i = 1:length(states)
        pressure_prof_3D(i + 2, 3:end-1) = states{i, 1}.pressure(2:end-1,1)';
        pressure_prof_3D(i + 2,2) = interp1(centroids(2:3,1),states{i, 1}.pressure(2:3,1)',0,'linear','extrap');
        pressure_prof_3D(i + 2,end) = interp1(centroids(end-2:end-1,1) - ...
            2 * centroids(1,1),states{i, 1}.pressure(end-2:end-1,1),sample_L,'linear','extrap');
    end    
end

model.history_match.Sw_profile_3D = sw_profile_3D;

sw_profile_2D = zeros(length(params.cumScheduleSteps) + 1, x_dims + 1);
sw_profile_2D(2:end,1) = params.cumScheduleSteps;
sw_profile_2D(1,2:end) = [0;centroids(2:x_dims-1,1) - 2*centroids(1,1);sample_L];
sw_profile_2D(2:end,2:end) = mean(reshape(sw_profile_3D(2:end,2:end), ...
    [length(params.cumScheduleSteps) x_dims length(centroids)/x_dims]), 3);
sw_profile_2D(2,2) = interp1(centroids(2:3,1),state0.s(2:3,1),0,'linear','extrap');
sw_profile_2D(2,end) = interp1(centroids(end-2:end-1,1) - ...
    2 * centroids(1,1),state0.s(end-2:end-1,1),sample_L,'linear','extrap');
for i = 1:length(states)
    sw_profile_2D(i + 2,2) = interp1(centroids(2:3,1),states{i, 1}.s(2:3,1)',0,'linear','extrap');
    sw_profile_2D(i + 2,end) = interp1(centroids(end-2:end-1,1) - ...
        2 * centroids(1,1),states{i, 1}.s(end-2:end-1,1),sample_L,'linear','extrap');
end
model.history_match.Sw_profile = sw_profile_2D;

model.history_match.pressure_prof = pressure_prof_3D;
endofsched_t = model.history_match.endofsched_t;
[~,idx] = intersect(params.cumScheduleSteps,endofsched_t,'stable');
model.history_match.endofsched_satprofile = [sw_profile_3D(1,:) ; sw_profile_3D(idx+1,:)];
model.history_match.endofsched_press_profile = [pressure_prof_3D(1,:) ; pressure_prof_3D(idx+1,:)];
model.history_match.calculatedp = model.dynamic.params.pDiff;
model.history_match.endofsched_dp = model.history_match.calculatedp(idx);
model.history_match.SwAvg = model.dynamic.params.SwAvg;
model.history_match.endofsched_swavg = model.dynamic.params.SwAvg(idx);
process     = model.experiment.process;
processName = lower(process.name);            
if(strcmpi(processName,'drainage'))
    model.history_match.Qp_net = model.dynamic.params.Qp_net(:,1);
    model.history_match.endofsched_Qp_net = model.history_match.Qp_net(idx);
elseif(strcmpi(processName,'imbibition'))
    model.history_match.Qp_net = model.dynamic.params.Qp_net(:,2);
    model.history_match.endofsched_Qp_net = model.history_match.Qp_net(idx);
end


