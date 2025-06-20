function model = Run3D(model, visualization)
%
% DESCRIPTION: runs the simulation for defined schedule
%
% SYNOPSIS:
%   model = Run3D(model, visualization)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - geometry: length and diameter of the core
%   - rock: rock properties like porosity and absolute permeability
%   - fluid: fluid properties like viscosities and densities
%   - process: type of the simulation - SS, USS or Centrifuge, drainage or
%   imbibition
%   - simulation: time stepping and grid cells information
%   - schedule: flooding or rotation schedule depending on the type of the
%   experiment
%   - observation: path to files containing the measured experimental data
%   used for comparison in plots or history matching
%   - experiment: saturation functions used for forward modeling
%   - plot: plotting options
%   - output: define the properties desired in the output
%   - history_match: gradient based and MCMC history matching options
%
% RETURNS:
%   model - adds the the following fields to the struct:
%   - twoPhaseOilWaterModel: struct with the porosity, permeability and satNum index
%   - state0: the initial state of the model
%   - state: all the states during simulation time steps
%   - dynamic: information saved from simulation time steps
%
% ----------------------------------
% (c) 2024
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
static = SetupStaticParams(model); % set static struct
model.static = static;
state0 = Initialize(model); % set initial condition
state0.facePressure = zeros(model.grid.G.faces.num, 1) * ...
    model.experiment.schedule.pout.value;
model.state0 = state0;
model.state = state0;
dynamic = SetupDynamicParams(model); % initialize dynamic struct 
model.dynamic = dynamic;
rampupsteps = model.simulation.rampupsteps.value;
%---------------------------------------  
scheduleTable = model.experiment.schedule.procedure;
dt = model.static.dt;
%---------------------------------------
G     = model.grid.G;
rock  = model.rock;
fluid = model.fluid;
model.twoPhaseOilWaterModel = TwoPhaseOilWaterModel(G, rock, fluid);
%---------------------------------------
for i = 1 : height(scheduleTable)
    params = dynamic.params;
    params.periodStart = [params.periodStart,scheduleTable{i,1}];
    params.periodEnd = [params.periodEnd,scheduleTable{i,2}];
    params.periodInterval = params.periodEnd - params.periodStart;
    scheduleRow = scheduleTable{i,:};
    model.bc = BoudaryConditions(model,scheduleRow);
    processType = model.experiment.process.type; 
    if(strcmpi(processType,'cent'))
        model.gravity = scheduleRow(3);
        model = SetupCentrifuge(model);
    else
        gravity off
    end
    twoPhaseOilWaterModel = model.twoPhaseOilWaterModel;
    %---------------------------------------
    periodInterval = params.periodInterval(end);
    n = periodInterval/dt;
    if isfield(model.simulation, 'load_from_sat_prof') && model.simulation.load_from_sat_prof
        try
            rampup = diff(model.experiment.observation.satProfile.table...
                (1:end,1));
%             rampup = [repmat(rampup(1)/4, 4, 1); rampup(2:end)];
        catch
            fprintf('No saturation profile data provided!')
        end
    else
        rampup = rampupTimesteps(periodInterval, periodInterval/n, ...
            rampupsteps);
    end
    schedule = simpleSchedule(rampup, 'bc', model.bc);
    params.scheduleSteps = rampup;
    params.cumScheduleSteps = [params.cumScheduleSteps;...
    cumsum(params.scheduleSteps) + params.cumScheduleSteps(end)];
    schedule.counter = i;
    % simulate the schedule
    %---------------------------------------
    verbose = false;
  [model.state, model.schedulereport] = Simulate(...
      model,twoPhaseOilWaterModel,schedule, verbose);
    % we use this to hide the time steps report from MRST
    % [~, model.state, model.schedulereport] = evalc...
    %     ('Simulate(model,twoPhaseOilWaterModel,schedule, verbose)');
    %---------------------------------------
    params.schedule = [params.schedule; schedule];
    params.fromIdx = params.toIdx;
    params.counter = params.counter + length(schedule.step.control);        
    params.toIdx   = params.counter;
    params.plot = false;
    dynamic.params = params;             
    dynamic.states = [dynamic.states;model.state]; 
    dynamic.schedulereport = [dynamic.schedulereport; model.schedulereport];
    model.dynamic  = dynamic;
    dynamic = UpdateDynamicParams_3D(model); % update dynamic parameters
    model.dynamic = dynamic;
    if visualization
        params.plotNo = params.plotNo + 1; 
        Visualize_3D(model); 
        frames(params.plotNo) = getframe(gcf); 
    end
    dynamic.params.plotNo = params.plotNo;
end  
if visualization
    model.dynamic.params.frames = frames;
    if not(isempty(model.output))
        if(model.output.include)
            SaveResults(model);
        end
    end
end
%---------------------------------------
% saving parameters for history match
%---------------------------------------
model = save_for_HM_3D(model,state0);
%---------------------------------------
% plotting diagnostics toolbar
figure
plotToolbar(model.grid.G, model.dynamic.states); colormap parula(20);
view(0,0)