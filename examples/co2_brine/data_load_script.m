[status,sheets] =xlsfinfo('06_ESTLC5H_unsatCO2brine_drainage_scan01_as_scan00_24x24x8_-180_corr.xlsx');
satprofiles_3d = cell(length(sheets),1);
for k=1:numel(sheets)
  satprofiles_3d{k} = readtable('06_ESTLC5H_unsatCO2brine_drainage_scan01_as_scan00_24x24x8_-180_corr.xlsx', 'Sheet', sheets{k}, 'ReadVariableNames',false);
end
