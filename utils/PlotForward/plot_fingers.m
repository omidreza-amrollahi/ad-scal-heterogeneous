function plot_fingers(model, time_index)
%
% DESCRIPTION: adds fluid model to the model struct
%
% SYNOPSIS:
%   plot_fingers(model, time_index)
%
% PARAMETERS:
%   model - struct output from the simulation
%   time_index - index of the simulation timestep to plot the finger snapshot
%
% RETURNS:
%   plot showing the viscous fingers
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

newSw = 1-model.history_match.Sw_profile_3D([1,time_index],2:end);
injtime = model.history_match.Sw_profile_3D(time_index,1);
reshaped_array = reshape(newSw(2,:), [model.grid.G.cartDims(1), numel(newSw(2,:))/model.grid.G.cartDims(1)])';

saturation_level = 1-0.95;
figure
imagesc(reshaped_array);
hold on;
contour(reshaped_array, [saturation_level, saturation_level], 'k', 'LineWidth', 1);

% find peaks in the contour line
C = contourc(reshaped_array, [saturation_level, saturation_level]);
peak_x = findpeaks(C(1,:) , 'MinPeakDistance', 10, 'MinPeakHeight',mean(C(1,:)));
peak_count = numel(peak_x);
[tf, idx] = ismember(peak_x, C(1,:));
peak_y = C(2,idx);

% plot markers at peak positions
filtered_columns = C(:, C(1,:) > 25);
plot(peak_x, peak_y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 3);
fprintf('mean location of front: %.2f.\n', mean(filtered_columns(1,1:end)) )
fprintf('Max-Min location of front: %.2f.\n', max(filtered_columns(1,1:end)) - min(filtered_columns(1,1:end)) )
% display peak count
% title(['Injection Time: ' num2str(injtime/3600) ' hours']);

title(['Number of peaks: ' num2str(peak_count) ' Injection Time: ' num2str(injtime/3600) ...
    ' hours, Max-Min location of front: ' num2str(max(filtered_columns(1,1:end)) - min(filtered_columns(1,1:end))) ...
    ' ,Mean location of front: ' num2str(mean(filtered_columns(1,1:end)))]);
disp(['Number of peaks: ' num2str(peak_count) ' Injection Time: ' num2str(injtime/3600) ' hours']);

set(gcf, 'Units', 'pixels');
width = 800;  % Specify the width of the figure (in pixels)
height = width / 2;  % Calculate the height based on the desired aspect ratio
set(gcf, 'Position', [100, 100, width, height]);
set(gcf, 'PaperPositionMode', 'auto');

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.XTickLabel = {};
ax.YTickLabel = {};

colormap parula(10);
clim([0, 0.5]);

% Save figure
% filename = sprintf([foldername,'/figure_%d.png'], time_index);
% saveas(gcf, filename, 'png');
% close(gcf);