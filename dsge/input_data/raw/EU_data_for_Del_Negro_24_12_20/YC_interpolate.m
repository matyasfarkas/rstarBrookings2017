% Load the data (Assume the data is in a table format with columns: Date, OIS3m, OIS6m, OIS1y, OIS2y)
data = readtable('OIS.csv'); % Replace with the actual filename

% Extract the monthly series
dates = data.date;
OIS3m = data.OIS3m;
OIS6m = data.OIS6m;
OIS1y = data.OIS1y;
OIS2y = data.OIS2y;

% Define maturities (in months) for interpolation
monthly_maturities = [3, 6, 12, 24]; % Original maturities in months
quarterly_maturities = [3, 6, 9, 12, 15, 18]; % Required maturities in months

% Preallocate for interpolated yields
interpolated_yields = zeros(length(dates), length(quarterly_maturities));

% Interpolate yields for each date
for i = 1:length(dates)
    if i == 25
        i  =i +1;
    end
    % Original data for the current date
    original_maturities = monthly_maturities;
    yields = [OIS3m(i), OIS6m(i), OIS1y(i), OIS2y(i)];
    
    % Interpolate for required maturities
    interpolated_yields(i, :) = interp1(original_maturities, yields, quarterly_maturities, 'cubic');
end

% Convert dates to quarters and aggregate using averages
quarterly_dates = dateshift(dates, 'start', 'quarter');
[unique_quarters, ~, quarter_idx] = unique(quarterly_dates);

% Preallocate for aggregated yields
aggregated_yields = zeros(length(unique_quarters), length(quarterly_maturities));

% Compute quarterly averages
for q = 1:length(unique_quarters)
    aggregated_yields(q, :) = mean(interpolated_yields(quarter_idx == q, :), 1);
end

% Combine into a table
result = array2timetable(aggregated_yields, ...
    'VariableNames', {'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6'},'RowTimes',unique_quarters);
% result.Quarters = unique_quarters;

% Display the final result
disp(result);

% Save the result to a CSV file
writetimetable(result, 'yield_curve_quarters.csv');
