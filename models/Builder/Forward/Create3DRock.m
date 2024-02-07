function model = Create3DRock(model)
% DESCRIPTION:
% Adds 3D gridding information to the model struct, particularly for
% creating a 3D grid based on the experiment geometry and simulation
% settings. This includes calculating grid cell dimensions and positions
% in the x, y, and z directions, and adjusting for boundary cells if
% specified. It also computes the geometry of the grid, identifies inlet
% and outlet cells, and optionally creates saturation number regions for
% simulation.

% SYNOPSIS:
%   model = Create3DGrid(model)

% PARAMETERS:
%   model - A struct containing the following fields:
%       - experiment: Contains saturation functions and geometry (length and diameter) for modeling.
%       - simulation: Includes time stepping, grid cells information, and optionally boundary cells.

% RETURNS:
%   model - The input struct updated with the following fields:
%       - grid: A struct with detailed information about the created 3D grid, including:
%           - area, dx, dy, dz: Area of grid cells and their dimensions.
%           - x, y, z: Arrays defining the positions of grid nodes.
%           - dxs, dys, dzs: Arrays defining the sizes of grid cells.
%           - G: The generated grid, after computing geometry and removing cells outside the specified domain.
%           - remove_mask: Indicates which cells were removed based on geometry.
%           - Inlet and outlet masks: Logical arrays identifying inlet and outlet cells for boundary conditions.
%           - satNum: (Optional) Array defining saturation number regions for simulation.

% COMMENTS:
% - The function handles both uniform and boundary-adjusted grid creation.
% - It is designed to support simulations with specific geometrical constraints, such as cylindrical domains.
% - The implementation assumes a cylindrical coordinate system with adjustments for numerical errors in identifying boundary cells.

% ----------------------------------
% (c) 2024
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
    experiment = model.experiment;
    simulation = model.simulation;
    grid       = model.grid;
    G          = grid.G;

    %%
    poro = experiment.rock.poro.value * ones(G.cells.num, 1);
%% sensitivity analysis for porosity
%     poro = bsxfun(@plus,poro,randn(numel(poro),1).*poro*0.5);
%     poro(poro < 0.001) = 0.001; poro(poro > 1) = 1;
%     poro(1:39*70) = experiment.rock.poro.value * 2;
%     poro(39*70+1:39*155) = experiment.rock.poro.value * 3;
%     figure
%     plotToolbar(G, poro);
%     title('Porosity')
    %%
    perm = experiment.rock.perm.value * ones(G.cells.num, 1);
%% sensitivity analysis for permeability
    % perm = perm + randn(numel(perm),1) .* perm * 0.04;
    rng(123);
    perm = perm + (rand(numel(perm),1)-0.5) .* perm * 0.01;
    % perm(perm < 0) = 1*milli*darcy;
%     perm(1:39*70) = experiment.rock.perm.value * 2;
%     perm(39*70+1:39*155) = experiment.rock.perm.value * 3;


    % figure
    % plotToolbar(G, perm/milli/darcy);
    % title('Permeability')

%% assigning random fields to porosity and permeability
%     poro = bsxfun(@plus,poro,randn(numel(poro),1).*poro*0.5);
%     poro(poro < 0.001) = 0.001; poro(poro > 1) = 1;
%     perm = bsxfun(@plus,perm,randn(numel(perm),1).*perm*0.3);
%     figure; plot(perm*1e14)
%     poro = use_porosity_profile(simulation, G, experiment);
%     figure; plot(G.cells.centroids(:,1), poro);
       
    rock = makeRock(G, perm, poro); 
    rock.pv = poreVolume(G, rock); 
    if (isfield(simulation,'bCells'))
        rock.regions = struct('saturation', grid.satNum);
    end
    model.rock = rock;



end