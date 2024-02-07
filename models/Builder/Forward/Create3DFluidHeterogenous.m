function model = Create3DFluidHeterogenous(model, f, sw_diff)

muW  = model.experiment.fluid.muW.value;
muO  = model.experiment.fluid.muNW.value;
rhoW = model.experiment.fluid.rhoW.value;
rhoO = model.experiment.fluid.rhoNW.value;
simulation = model.simulation;

fluid = initSimpleADIFluid('phases', 'WO', ...
                           'mu'    , [muW, muO]);
fluid.isIncomp = true;
fluid.rhoWS = rhoW;
fluid.rhoOS = rhoO;  

if (isfield(simulation,'bCells'))   
    fluid.krW = cell(1, length(f) + 1);
    [fluid.krW{1, :}] = deal(@(sw) interpTable(model.satfun.sw_kr, model.satfun.krw, sw));
    fluid.krW{1, 1} = @(sw) sw;

    fluid.krO = cell(1, length(f) + 1);
    [fluid.krO{1, :}] = deal(@(so) interpTable(model.satfun.sw_kr, model.satfun.kro, 1-so));
    fluid.krO{1, 1} = @(so) so;

    fluid.pcOW = cell(1, length(f) + 1);
    fluid.pcOW{1, 1} = @(s) 0;
    for i = 2: length(f) + 1
        fluid.pcOW{1, i} = @(sw) interpTable(model.satfun.sw_pc - sw_diff(i-1), model.satfun.pc./f(i-1), sw);
    end

end    
model.fluid = fluid;
