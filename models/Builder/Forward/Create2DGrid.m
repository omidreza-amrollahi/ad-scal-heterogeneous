function model = Create2DGrid(model)
%
% DESCRIPTION: adds griding information to the model struct
%
% SYNOPSIS:
%   model = CreateGrid(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling
%   - simulation: time stepping and grid cells information
%
% RETURNS:
%   model - adds the the following fields to the struct:
%   - grid: struct with the information about gridding
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
% 2-D grid in x-direction
experiment = model.experiment;
simulation = model.simulation;
length = experiment.geometry.length.value;
diameter = experiment.geometry.diameter.value;
nCells = simulation.nCells.value;
dx = length / nCells;
area = pi * diameter ^ 2 / 4;
model.grid.area = area;
dy = sqrt(area); dz = dy;
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
    y = [0;dy]; model.grid.y = y;
    z = linspace(0, dz, number_of_y_cells); model.grid.z = z;
    dxs = [dxLeft,dx * ones(1,nCells),dxRight]; model.grid.dxs = dxs;
    dys = dy * ones(nCells+2,1); model.grid.dys = dys;
    dzs = dz * ones(nCells+2,1); model.grid.dxz = dzs;        
else
    x = (0:dx:nCells*dx)'; model.grid.x = x;
    y = [0;dy]; model.grid.y = y;
    z = linspace(0, dz, number_of_y_cells); model.grid.z = z;
    dxs = dx * ones(nCells,1); model.grid.dxs = dxs;
    dys = dy * ones(nCells,1); model.grid.dys = dys;
    dzs = dz * ones(nCells,1); model.grid.dzs = dzs;
end

% create MRST grid
G = tensorGrid(x, y, z);
G = computeGeometry(G);
G = gridAddHelpers(G);
grid.G = G;  

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

% Display the grid
%---------------------------------------
%     figTitle   = 'Grid';
%     xLabel     = 'x';
%     yLabel     = 'y'; 
%     zLabel     = 'z';
%     figStyle   = 'docked';
%     figTag     = figTitle;
%     fig = figure('Name',        figTitle, ...
%                  'Tag',         figTag,   ...
%                  'NumberTitle', 'off',     ...
%                  'WindowStyle', figStyle );
%     figure(fig);
%     xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
%     xLim = [min(x) max(x)]; yLim = [min(y) max(y)]; zLim = [min(z) max(z)];
%     xlim(xLim); ylim(yLim); zlim(zLim);
%     plotGrid(G,'FaceColor',[.7 .7 1]); view(3); axis equal;
%     title('Grid Representation');
%---------------------------------------

model.grid = grid;