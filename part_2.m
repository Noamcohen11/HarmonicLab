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

%resutls = {path, mass in kg}
lab_results = {
    'part_2\0.5.csv'            , 0.5;
    'part_2\3.68.csv'           , 3.68;
    'part_2\3.72.csv'           , 3.72;
    'part_2\1.5.csv'            , 1.5;
    'part_2\1.csv'              , 1;
    'part_2\2.csv'              , 2;
    'part_2\3.5.csv'            , 3.5;
    'part_2\3.7.csv'            , 3.7;
    'part_2\3.8.csv'            , 3.8;
    'part_2\3.9.csv'            , 3.9;
    'part_2\3.65.csv'           , 3.65;
    'part_2\3.73.csv'           , 3.73;
    'part_2\3.74.csv'           , 3.74;
    'part_2\3.75.csv'           , 3.75;
    'part_2\3.85.csv'           , 3.85;
    'part_2\3.87.csv'           , 3.87;
    'part_2\3.csv'              , 3;
    'part_2\4.2.csv'            , 4.2;
    'part_2\4.5.csv'            , 4.5;
    'part_2\4.csv'              , 4;
    'part_2\6.csv'              , 6;
    };

%% code:
% Fit each result:
aplitude = zeros(1,size(lab_results,1));
herz = zeros(1,size(lab_results,1));
for i = 1:size(lab_results,1)
    results = readtable(string(lab_results(i,1)));
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
        amplitude_fit = DampedAmplitudeFit(locs,pks);
        peak_error = zeros(1,length(locs)) + distance_error_range;
        time_error = zeros(1,length(locs));
        aplitude(i) = abs(pks(length(pks)));
    end
end

%% Plot for all mass
figure
hold on
herz = cell2mat(lab_results(:,2))';
herz_error = zeros(1,length(herz));
amplitude_error = zeros(1,length(aplitude)) + distance_error_range;
errorbar(herz , aplitude, herz_error, herz_error, amplitude_error, amplitude_error, 'color','magenta','LineStyle','none', 'LineWidth', 2)
dumbell_fit = Dumbellfit(herz' ,aplitude');
plot(dumbell_fit)
xlabel('Frequency(hz)')
ylabel('amplitude(m)')
legend('peaks', 'Fitted Curve')
hold off

function f = DampedAmplitudeFit(x, y)
    fitt = fittype('a*exp(-b*x) -c*x','coefficients', {'a', 'b', 'c'});
    f = fit(x,y,fitt);
end

function f = Dumbellfit(x, y)
    fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', [0.008463, 14.34, 0.09421]);
    fitt = fittype('a/(sqrt((b-x^2)^2 + (c*x)^2))','coefficients', {'a', 'b', 'c'}, 'options', fo);
    f = fit(x,y,fitt);
end


