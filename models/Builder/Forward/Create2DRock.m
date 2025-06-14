function model = Create2DRock(model)
%
% DESCRIPTION: adds porosity and pearmeability information to the model
%
% SYNOPSIS:
%   model = CreateRock(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling
%   - simulation: time stepping and grid cells information
%   - grid: struct with the information about gridding
%
% RETURNS:
%   model - adds the the following fields to the struct:
%   - rock: struct with the porosity, permeability and satNum index
%
% ----------------------------------
% (c) 2023
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
experiment = model.experiment;
simulation = model.simulation;
grid       = model.grid;
G          = grid.G;
poro = experiment.rock.poro.value * ones(G.cells.num, 1);
perm = experiment.rock.perm.value * ones(G.cells.num, 1);
perm_het = experiment.rock.perm_heterogeneity.include;
perm_het_value = experiment.rock.perm_heterogeneity.value;

if perm_het
    rng(123)
    perm = perm + (rand(numel(perm),1)-0.5) .* perm * ...
    perm_het_value;
end

rock = makeRock(G, perm, poro); 
rock.pv = poreVolume(G, rock); 
if (isfield(simulation,'bCells'))
    rock.regions = struct('saturation', grid.satNum);
end
model.rock = rock;