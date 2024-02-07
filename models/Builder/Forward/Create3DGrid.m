function model = Create3DGrid(model)
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
    % 3-D grid in x-direction
    experiment = model.experiment;
    simulation = model.simulation;
    length = experiment.geometry.length.value;
    diameter = experiment.geometry.diameter.value;
    nCells = simulation.nCells.value;
    dx = length / nCells;
    area = pi * diameter ^ 2 / 4;
    model.grid.area = area;
    % note this part is different from 1D version
    dy = diameter; dz = dy;
    model.grid.dx = dx; model.grid.dy = dy; model.grid.dz = dz;
    % number_of_y_cells is the number of cells in the y-z direction
    % +1 is necessary for proper functionality
    number_of_y_cells = model.simulation.nCells_y.value + 1;
    if (isfield(simulation,'bCells'))        
        bCells = simulation.bCells.value;
        x = zeros(nCells+3,1);
        dxLeft  = bCells * dx;
        dxRight = bCells * dx;
        x(1) = 0;
        x(2) = dxLeft; 
        x(3:end-1) = dxLeft + (dx:dx:nCells*dx)';
        x(end) = x(end-1) + dxRight;
        model.grid.x = x;
        y = linspace(0, dy, number_of_y_cells); model.grid.y = y;
        z = linspace(0, dz, number_of_y_cells); model.grid.z = z;
        dxs = [dxLeft,dx * ones(1,nCells),dxRight]; model.grid.dxs = dxs;
        dys = dy * ones(nCells+2,1); model.grid.dys = dys;
        dzs = dz * ones(nCells+2,1); model.grid.dxz = dzs;        
    else
        x = (0:dx:nCells*dx)'; model.grid.x = x;
        y = linspace(0, dy, number_of_y_cells); model.grid.y = y;
        z = linspace(0, dz, number_of_y_cells); model.grid.z = z;
        dxs = dx * ones(nCells,1); model.grid.dxs = dxs;
        dys = dy * ones(nCells,1); model.grid.dys = dys;
        dzs = dz * ones(nCells,1); model.grid.dzs = dzs;
    end
   
    % create MRST grid
%     G = cartGrid([nCells+2,1,1], [max(x),max(y),max(z)]*meter);    
    G = tensorGrid(x, y, z);
    G = computeGeometry(G);
    G = gridAddHelpers(G);
    
    radius = diameter/2;
    euclidean_distances = sqrt((G.cells.centroids(:,2)-radius).^2+(G.cells.centroids(:,3)-radius).^2);
    remove_mask = euclidean_distances > radius;
    G = removeCells(G, remove_mask);
    grid.G = G;
    grid.remove_mask = remove_mask;

    % find index of inlet face
    % this way is for numerical errors
    eps = 1e-8;
    inlet_mask = abs(G.cells.centroids(:,1) - G.cells.centroids(1,1)) < eps;
    outlet_mask = abs(G.cells.centroids(:,1) - G.cells.centroids(G.cartDims(1,1), 1)) < eps;
    inlet_frontcells_mask = abs(G.cells.centroids(:,1) - G.cells.centroids(2,1)) < eps;
    outlet_backcells_mask = abs(G.cells.centroids(:,1) - G.cells.centroids(G.cartDims(1,1) - 1, 1)) < eps;

    grid.G.inlet_mask = inlet_mask;
    grid.G.outlet_mask = outlet_mask;
    grid.G.inlet_frontcells_mask = inlet_frontcells_mask;
    grid.G.outlet_backcells_mask = outlet_backcells_mask;
    
    % create satnum regions
    if (isfield(simulation,'bCells'))
        satNum       = ones(grid.G.cells.num, 1) * 2;
        satNum(inlet_mask)    = 1; satNum(outlet_mask) = 1;
        grid.satNum  = satNum;     
    end

%     Display the grid
% ------------------------------------------------------------------------
    % figTitle   = 'Grid';
    % xLabel     = 'x';
    % yLabel     = 'y'; 
    % zLabel     = 'z';
    % figStyle   = 'docked';
    % figTag     = figTitle;
    % fig = figure('Name',        figTitle, ...
    %              'Tag',         figTag,   ...
    %              'NumberTitle', 'off',     ...
    %              'WindowStyle', figStyle );
    % figure(fig);
    % xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
    % xLim = [min(x) max(x)]; yLim = [min(y) max(y)]; zLim = [min(z) max(z)];
    % xlim(xLim); ylim(yLim); zlim(zLim);
    % plotGrid(G,'FaceColor',[.7 .7 1]); view(3); axis equal;
    % title('Grid Representation');
% ------------------------------------------------------------------------
%     figure
%     plotToolbar(G, [inlet_mask, outlet_mask, satNum]);

    model.grid = grid;
end