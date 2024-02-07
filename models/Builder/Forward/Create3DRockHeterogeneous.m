function model = Create3DRockHeterogeneous(model, porosity_profile, f_factor)
    experiment = model.experiment;
    simulation = model.simulation;
    grid       = model.grid;
    G          = grid.G;
    coreLength = model.experiment.geometry.length.value;
    poro_array = experiment.rock.poro.value * ones(G.cells.num, 1);
    perm_array = experiment.rock.perm.value * ones(G.cells.num, 1);
  
    if (isfield(simulation,'bCells'))
%         if isfield(experiment.rock.poro, 'porosity_profile')
%             porosity_profile_obs = experiment.rock.poro.porosity_profile;
%             length_obs = porosity_profile_obs(:,1); % care for the units
%             porosity_obs = porosity_profile_obs(:,2);
%             x = [0;G.cells.centroids(2:end-1,1) - 2 * G.cells.centroids(1,1);coreLength];
%             poro = interp1(length_obs / 100, porosity_obs, x); 
%         end
        % figure; plot(x,poro); xlabel("Distance (m)"); ylabel("Porosity")
%--
        inner_mask = not(model.grid.G.inlet_mask) & not(model.grid.G.outlet_mask);
        PoroInputProfile = porosity_profile.porosity_profile_vertical;
        PoroInputProfile(PoroInputProfile==0) = [];
        poro_array(inner_mask) = PoroInputProfile;
        perm_array(inner_mask) = perm_array(inner_mask) ./ f_factor .^2 .* poro_array(inner_mask) ./ experiment.rock.poro.value;
%--
        figure; plot(perm_array); harmmean(perm_array);
    end

    rock = makeRock(G, perm_array, poro_array); 
    rock.pv = poreVolume(G, rock); 
    if (isfield(simulation,'bCells'))
        rock.regions = struct('saturation', grid.satNum);
    end
    model.rock = rock;
end