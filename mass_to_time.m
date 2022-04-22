%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 20/04/2022   %
%   Lab - experiment 2      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mass to time summary:
% Take damped harmonic data and fit it as:
% a*exp(b*x)*(cos(2*pi*x/c + 2*pi/d))+e
% Prints goodness of fit and plot

%% Parameters:
plot_every_mass = 1;
distance_error_range = 0.0172;
lab_results = {
    'csv_files\5 G part 2.csv'     , 0.005;
    'csv_files\10 G part 2.csv'    , 0.01;
    'csv_files\14.6 G part 2.csv'  , 0.0146;
    'csv_files\20 G part 1.csv'    , 0.02;
    %'csv_files\84 G part 2.csv'    , 0.0847;
    };

%% code:
% Fit each result:

CycleTime = zeros(1,size(lab_results,1));

for i = 1:size(lab_results,1)

    %% Grab lab lab_results
    results = readtable(string(lab_results(i,1)));
    y = results{:,2};
    x = results{:,1};

    %% Remove unwanted data

    % Remove data captured before the experiment begins:
    init_wave = find(y == min(y));
    x = x-x(init_wave);
    for j = 1:init_wave-1
        x(1) = [];
        y(1) = [];
    end
    
    %Remove data captured after the experiment ends:
    act_zero = abs(y(length(y)))+distance_error_range;
    for j = 60:length(x)
        if (abs(y(j))<act_zero) && (abs(y(j-60))<act_zero) 
            x(j:length(x)) = [];
            y(j:length(y)) = [];
        break
        end
    end
    
    final_fit = DampedHarmonic_fit(x,y);
    % Find peaks:
    [pks,locs] = findpeaks(y,x);
    amplitude_fit = DampedAmplitudeFit(locs,pks);
    peak_error = zeros(1,length(locs)) + distance_error_range;
    time_error = zeros(1,length(locs));
    CycleTime(i) = 2*pi/mean(diff(locs(2:10)));
    %fit_values = coeffvalues(final_fit);
    %CycleTime(i) = fit_values(3);
    
    %% Plot per mass:
    if plot_every_mass
        figure
        % Plot actual points:
        plot(x,y, '.')
        hold on
        % Plot peaks
        plot(final_fit,locs, pks, 'o')
        hold off
        grid
        box on
        xlabel('Time(S)')
        ylabel('Distance(M)')
        legend('Original Data', 'peaks', 'Fitted Curve')
        
        %Plot peaks per time
        figure
        hold on
        errorbar(locs, pks, peak_error, peak_error, time_error, time_error, 'color','magenta','LineStyle','none', 'LineWidth', 2)
        plot(amplitude_fit)
        grid
        box on
        xlabel('peaks time(S)')
        ylabel('peaks(M)')
        legend('peaks', 'Fitted Curve')
        hold off
    end
end


%% Plot for all mass
figure
mass = cell2mat(lab_results(:,2))';
plot(cell2mat(lab_results(:,2))',2*pi./CycleTime, '.')
hold on
xlabel('Mass(KG)')
ylabel('omega(Rad/S)')
legend('peaks', 'Fitted Curve')
hold off

function f = DampedHarmonic_fit(x, y)
    %% Get fit parametes:
    
    y = detrend(y);                                                                      % Remove Linear Trend
    yu = max(y);
    yl = min(y);
    yr = (yu-yl);                                                                        % Range of  y’
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                                 % Returns Approximate Zero-Crossing Indices Of Argument Vector
    zt = x(zci(y));
    per = 2*mean(diff(zt));                                                              % Estimate period
    ym = mean(y);                                                                        % Estimate offset
    
    init_fit = @(b,x)  b(1) .* exp(b(2).*x) .* (cos(2*pi*x./b(3) + 2*pi/b(4))) + b(5);   % Objective Function to fit
    fcn = @(b) norm(init_fit(b,x) - y);                                                  % Least-Squares cost function
    [s,] = fminsearch(fcn, [yr; -10;  per;  -1;  ym]);                                   % Minimise Least-Squares
    fit_params = s';
    %% Fit:
    fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', fit_params);         % Use the parameters gathered as starting points.
    fitt = fittype('a.*exp(b*x).*(cos(2*pi.*x/c + 2.*pi./d))+e','coefficients', {'a', 'b', 'c', 'd', 'e'}, 'options', fo);
    f = fit(x,y,fitt);
end

function f = DampedAmplitudeFit(x, y)
    fitt = fittype('a*exp(-b*x)+c','coefficients', {'a', 'b', 'c'});
    f = fit(x,y,fitt);
end

