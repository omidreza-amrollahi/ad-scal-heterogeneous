
%% load the satprofiles_3d
load('satprofiles_3d.mat')
load('time_array.mat')

%% process and prepare the satprofiles_3d
% +1 because satprofiles_3d are comming from python and there index starts from 0
max_x = max(satprofiles_3d{1}.Var1) + 1;
min_x = min(satprofiles_3d{1}.Var1) + 1;
max_y = max(satprofiles_3d{1}.Var2) + 1;
min_y = min(satprofiles_3d{1}.Var2) + 1;
max_z = width(satprofiles_3d{1}) - 2;

offset = -0.0085; 

% k is for time
for k = 1:length(satprofiles_3d)
    x_vector = satprofiles_3d{k}.Var1(3:end) + 1;
    y_vector = satprofiles_3d{k}.Var2(3:end) + 1;
    z_matrix_2d = satprofiles_3d{k}{3:end, 3:end};
    size_z_vector = size(z_matrix_2d);
    
    z_matrix_3d = zeros(size_z_vector(2), max_x , max_y) - offset;
    % j is for each slice
    for j = 1:size_z_vector(2)
        % i is for y,x axis
        for i=1:length(y_vector)
            z_matrix_3d(j, x_vector(i) , y_vector(i)) = z_matrix_2d(i, j);
        end
    end

%     z_matrix_2d = imresize(z_matrix_2d , 0.125, "bilinear");
%     z_matrix_3d = imresize3(z_matrix_3d , [37 17 17], "triangle");
    mean_1d{k} = mean(z_matrix_2d, 1)' + offset;
    % because in the incomming satprofiles_3d, the core length is in the z axis but
    % the simulator is coded in a way that the core length is in the x axis
    z_matrix_3d_saved{k} = z_matrix_3d + offset;
    z_matrix_vertical{k} = z_matrix_3d(:) + offset;
    mean_new{k} = mean(mean(z_matrix_3d_saved{k},2),3);

end

%% plotting using MRST tools
CoreLength = 15; %cm
CoreDiameter = 7.5; %cm
size_z_matrix_3d = size(z_matrix_3d);
a = linspace(0, CoreDiameter, size_z_matrix_3d(2) + 1);
b = linspace(0, CoreDiameter, size_z_matrix_3d(3) + 1);
c = linspace(0, CoreLength, size_z_matrix_3d(1) + 1);
G = tensorGrid(c, b, a);
G = computeGeometry(G);
figure
h1 = plotToolbar(G, struct('z_matrix_vertical', z_matrix_vertical));
xlabel('Core Length'); ylabel('y axis'); zlabel('z axis');
grid on
colormap hot
colorbar
xlim([0 CoreLength]); ylim([0 CoreDiameter]); zlim([0 CoreDiameter]);
title("CO2 Saturation")

%% validated with previous 1d observations
SatProfile1DLocation = model.experiment.observation.satProfile.table(1 ,2:end);
SatProfile1D= model.experiment.observation.satProfile.table(2:end,2:end);
SatProfile1DTime = model.experiment.observation.satProfile.table(2:end, 1);
MeasuredSlicesN = length(SatProfile1D);
time_array_s = time_array*60; % Assume time_array is in minutes
mutual_idx_t_obs = ismember(SatProfile1DTime,time_array_s, 'row');
SatProfile1D = SatProfile1D(mutual_idx_t_obs, :);
mutual_idx_t = ismember(time_array_s,SatProfile1DTime, 'row');
SliceAverageOf3DSatProf = mean_1d(mutual_idx_t);
figure
tiledlayout("flow")
ax1 = nexttile([1 2]);
time_index = [4, 22];
p1 = plot(ax1, SatProfile1DLocation, SatProfile1D(time_index, :));
hold(ax1, "on");
SatProfile1DInterped = ones(length(time_index), size_z_matrix_3d(1));
SatProfile3DtoPlot = ones(length(time_index), size_z_matrix_3d(1));
for i = 1:length(time_index)
    SatProfile1DInterped(i,:) = interp1(linspace(0,15,MeasuredSlicesN), SatProfile1D(time_index(i),:), linspace(0,15,size_z_matrix_3d(1)));
    SatProfile3DtoPlot(i,:) = (1 - SliceAverageOf3DSatProf{1, time_index(i)})';
end
SatProfile3DtoPlot = SatProfile3DtoPlot'; SatProfile1DInterped = SatProfile1DInterped'; 
% fun = @(coeffs) objective_function(coeffs, SatProfile3DtoPlot, SatProfile1DInterped);
% x0 = ones(1,100);
% options = optimoptions('fmincon','Display','iter');
% x = fmincon(fun, x0, [],[],[],[],[],[],[],options);
% SatProfile3DtoPlotFitted = polynomial(x, SatProfile3DtoPlot);
SliceLocationInterped = interp1(linspace(0,CoreLength,MeasuredSlicesN), SatProfile1DLocation, linspace(0,CoreLength,size_z_matrix_3d(1)));
ylim([0.4 1])
p2 = plot(ax1, SliceLocationInterped, SatProfile3DtoPlot, 'DisplayName', 'satprofiles_3d');
for i=1:length(p1)
    p2(i).Color = p1(i).Color;
end
nexttile;
difference = mean(SatProfile3DtoPlot - SatProfile1DInterped,2);
s = scatter(SliceLocationInterped, difference);
diff_mean = mean(difference);
diff_std = std(difference);
l1 = yline(diff_mean, 'k--', 'DisplayName', 'Mean');
l2 = yline(diff_mean+diff_std, 'r--', 'DisplayName', 'Std');
yline(diff_mean-diff_std, 'r--')
title(sprintf("Difference between two satprofiles_3d with mean %.4f",diff_mean))
legend([l1, l2])
ax3 = nexttile;
histfit(ax3, difference, 10, 'kernel')
set(ax3,'view',[90 -90])



%% apply the fitting results back to satprofiles_3d
z_matrix_3d_saved_new = cell(length(satprofiles_3d), 1);
z_matrix_vertical_new = cell(length(satprofiles_3d), 1);
for k = 1:length(satprofiles_3d)
    z_matrix_3d_saved_new{k} = polynomial(x, z_matrix_3d_saved{k});
    z_matrix_vertical_new{k} = polynomial(x, z_matrix_vertical{k});
end

%%
figure
h1 = plotToolbar(G, struct('z_matrix_vertical', z_matrix_vertical_new));
xlabel('Core Length'); ylabel('y axis'); zlabel('z axis');
grid on
colormap hot
colorbar
xlim([0 15]); ylim([0 7.5]); zlim([0 7.5]);
title("Decane Saturation")
%%
new_mean = cell(length(satprofiles_3d), 1);
for k = 1:length(satprofiles_3d)
    new_mean{k} =  mean(mean(z_matrix_3d_saved_new{k},2),3);
end

%% validated with previous 1d observations
SatProfile1DLocation = model.experiment.observation.satProfile.table(1 ,2:end);
SatProfile1D = model.experiment.observation.satProfile.table(2:end,2:end);
SatProfile1DTime = model.experiment.observation.satProfile.table(2:end, 1);
time_array_s = time_array*60;
mutual_idx_t_obs = ismember(SatProfile1DTime,time_array_s, 'row');
SatProfile1D = SatProfile1D(mutual_idx_t_obs, :);
mutual_idx_t = ismember(time_array_s,SatProfile1DTime, 'row');
SliceAverageOf3DSatProf = new_mean(mutual_idx_t);
figure
tiledlayout("flow")
ax1 = nexttile([1 2]);
p1 = plot(ax1, SatProfile1DLocation, SatProfile1D);
hold(ax1, "on");
time_index = 1:30;
SatProfile1DInterped = ones(length(time_index), 37);
sw_mean_1d = ones(length(time_index), 37);
for i = 1:length(time_index)
    SatProfile1DInterped(i,:) = interp1(linspace(0,15,296), SatProfile1D(time_index(i),:), linspace(0,15,37));
    sw_mean_1d(i,:) = (1 - SliceAverageOf3DSatProf{time_index(i),1})';
end
sw_mean_1d = sw_mean_1d'; SatProfile1DInterped = SatProfile1DInterped'; 
fun = @(coeffs) objective_function(coeffs, sw_mean_1d, SatProfile1DInterped);
x0 = ones(1,100);
options = optimoptions('fmincon','Display','iter');
% x = fmincon(fun, x0, [],[],[],[],[],[],[],options);
% satprofiles_3d_satProf = polynomial(x, sw_mean_1d);
SatProfile3DtoPlot = sw_mean_1d;
SliceLocationInterped = interp1(linspace(0,15,296), SatProfile1DLocation, linspace(0,15,37));
ylim([0.5 1])
p2 = plot(ax1, SliceLocationInterped, SatProfile3DtoPlot, 'DisplayName', 'satprofiles_3d');
for i=1:length(p1)
    p2(i).Color = p1(i).Color;
end
nexttile;
difference = mean(SatProfile3DtoPlot - SatProfile1DInterped,2);
s = scatter(SliceLocationInterped, difference);
diff_mean = mean(difference);
diff_std = std(difference);
l1 = yline(diff_mean, 'k--', 'DisplayName', 'Mean');
l2 = yline(diff_mean+diff_std, 'r--', 'DisplayName', 'Std');
yline(diff_mean-diff_std, 'r--')
title(sprintf("Difference between two satprofiles_3d with mean %.4f",diff_mean))
legend([l1, l2])
ax3 = nexttile;
histfit(ax3, difference, 10, 'kernel')
set(ax3,'view',[90 -90])

%% useful functions

function error = objective_function(coeffs, model, satprofiles_3d)
    diff = polynomial(coeffs, model) - satprofiles_3d;
    error = rms(diff, 'all');
end

function result = polynomial(coeffs, satprofiles_3d)
    result = zeros(size(satprofiles_3d));
    for i = length(coeffs):-1:2
        result = result + coeffs(i) * satprofiles_3d.^(i-1);
    end
end
