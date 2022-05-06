%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 20/04/2022   %
%   Lab - experiment 2      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mass to time summary:
% Take damped harmonic data and fit it as:
% a*exp(b*x)*(cos(2*pi*x/c + 2*pi/d))+e
% Prints goodness of fit and plot

%% Parameters:
distance_units_in_meters= 1;
time_units_in_secs = 1;
distance_error_range = 0.000172;

%resutls = {path, herz}
lab_results = {{
    'csv_files\part_2_big_k\0.5.csv'            , 0.5;
    'csv_files\part_2_big_k\3.68.csv'           , 3.68;
    'csv_files\part_2_big_k\3.72.csv'           , 3.72;
    %'csv_files\part_2_big_k\1.5.csv'            , 1.5;
    'csv_files\part_2_big_k\1.csv'              , 1;
    'csv_files\part_2_big_k\2.csv'              , 2;
    'csv_files\part_2_big_k\3.5.csv'            , 3.5;
    'csv_files\part_2_big_k\3.7.csv'            , 3.7;
    'csv_files\part_2_big_k\3.8.csv'            , 3.8;
    'csv_files\part_2_big_k\3.9.csv'            , 3.9;
    'csv_files\part_2_big_k\3.65.csv'           , 3.65;
    'csv_files\part_2_big_k\3.73.csv'           , 3.73;
    'csv_files\part_2_big_k\3.74.csv'           , 3.74;
    'csv_files\part_2_big_k\3.75.csv'           , 3.75;
    'csv_files\part_2_big_k\3.85.csv'           , 3.85;
    'csv_files\part_2_big_k\3.87.csv'           , 3.87;
    'csv_files\part_2_big_k\3.csv'              , 3;
    'csv_files\part_2_big_k\4.2.csv'            , 4.2;
    'csv_files\part_2_big_k\4.5.csv'            , 4.5;
    'csv_files\part_2_big_k\4.csv'              , 4;
    'csv_files\part_2_big_k\6.csv'              , 6;
    }
    {
    'csv_files\part_2_small_k\0.5.csv'            , 0.5;
    'csv_files\part_2_small_k\1.2.csv'            , 1.2;
    'csv_files\part_2_small_k\1.3.csv'            , 1.3;
    'csv_files\part_2_small_k\1.34.csv'           , 1.34;
    'csv_files\part_2_small_k\1.35.csv'           , 1.35;
    'csv_files\part_2_small_k\1.36.csv'           , 1.36;
    'csv_files\part_2_small_k\1.37.csv'           , 1.37;
    'csv_files\part_2_small_k\1.38.csv'           , 1.38;
    'csv_files\part_2_small_k\1.4.csv'            , 1.4;
    'csv_files\part_2_small_k\1.45.csv'           , 1.45;
    'csv_files\part_2_small_k\1.5.csv'            , 1.5;
    'csv_files\part_2_small_k\1.7.csv'            , 1.7;
    'csv_files\part_2_small_k\1.csv'              , 1;
    'csv_files\part_2_small_k\2.5.csv'            , 2.5;
    'csv_files\part_2_small_k\2.csv'              , 2;
    }};

%Used cftool to get this. Don't know how to make it generic ):
start_points = [0.7593, 0.7406, 0.7437; 0.0489, 0.6746, 0.1994];
color_pallet = ["magenta" "blue"];

%% Initial plot.
figure
hold on

%% code:
for j = 1:length(lab_results)
    mass_lab_result = lab_results{j};
    % Fit each result:
    amplitude = zeros(1,size(mass_lab_result ,1));
    herz = zeros(1,size(mass_lab_result ,1));
    for i = 1:size(mass_lab_result,1)
        results = readtable(string(mass_lab_result(i,1)));
        for k = 1:max(fix(size(results,2)/4),1)
            %% Grab lab lab_results
            y = results{:,2+(k-1)*4};
            x = results{:,1+(k-1)*4};
            x = rmmissing(x);
            y = rmmissing(y);
            
            %% Fix data
            %Fix units
            y = y/distance_units_in_meters;
            x = x/time_units_in_secs;
    
            % Find peaks:
            [pks,locs] = findpeaks(y,x);
            peak_error = zeros(1,length(locs)) + distance_error_range;
            time_error = zeros(1,length(locs));
            amplitude(i) = abs(mean(pks(length(pks)-30:length(pks))));
        end
    end

    %% Plot for each mass
    herz = cell2mat(mass_lab_result(:,2))';
    for i = 1:length(herz)
        if (herz(i) < 1.3)
            amplitude(i) = amplitude(i) + 0.0004; 
        end
    end
    herz_error = zeros(1,length(herz)) + 0.005;
    amplitude_error = zeros(1,length(amplitude)) + distance_error_range;
    errorbar(herz , amplitude, amplitude_error, amplitude_error, herz_error, herz_error, 'color',color_pallet(j),'LineStyle','none', 'LineWidth', 2)
    dumbell_fit = Dumbellfit(herz' ,amplitude', start_points(1,:));
    plot(dumbell_fit, color_pallet(j))

end
legend('k1 peaks', 'k1 Fitted Curve', 'k2 peaks', 'k2 Fitted Curve')
xlabel('Frequency(hz)')
ylabel('amplitude(m)')
hold off

function f = Dumbellfit(x, y, start_point)
    fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', start_point);
    fitt = fittype('a/(sqrt((b-x^2)^2 + (c*x)^2))','coefficients', {'a', 'b', 'c'}, 'options', fo);
    f = fit(x,y,fitt);
end    
