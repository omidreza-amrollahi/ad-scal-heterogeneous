function model = Create3DGridHeterogeneous(model)
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
    
    
    % find index of inlet face
    inlet_mask = G.cells.centroids(:,1) == G.cells.centroids(1,1);
    outlet_mask = G.cells.centroids(:,1) == G.cells.centroids(G.cartDims(1,1), 1);
    inlet_frontcells_mask = G.cells.centroids(:,1) == G.cells.centroids(2,1);
    outlet_backcells_mask = G.cells.centroids(:,1) == G.cells.centroids(G.cartDims(1,1) - 1, 1);

    grid.G.inlet_mask = inlet_mask;
    grid.G.outlet_mask = outlet_mask;
    grid.G.inlet_frontcells_mask = inlet_frontcells_mask;
    grid.G.outlet_backcells_mask = outlet_backcells_mask;

    % create satnum regions
    if (isfield(simulation,'bCells'))
        inner_mask = not(model.grid.G.inlet_mask) & not(model.grid.G.outlet_mask);
        satNum = (1:grid.G.cells.num)';
        satNum(or(inlet_mask, outlet_mask)) = 1;
        satNum(inner_mask) = 2:sum(inner_mask)+1;
        grid.satNum  = satNum;     
    else
        % not implemented
    end
    model.grid = grid;

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
% ------------------------------------------------------------------------
%     figure
%     plotToolbar(G, [inlet_mask, outlet_mask, grid.satNum]);
% ------------------------------------------------------------------------
end