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
distance_units_in_meters= 1;
time_units_in_secs = 1;
distance_error_range = 0.000172;
image_save_path = 'G:\My Drive\results\';
add_linear_damp = 0;

%resutls = {path, mass in kg}
lab_results = {
       'csv_files\5 G part 2.csv'     , 0.0048;
       'csv_files\10 G part 2.csv'    , 0.0098;
       'csv_files\14.6 G part 2.csv'  , 0.0146;
       'csv_files\20 G part 1.csv'    , 0.0199;
       'csv_files\25G.csv'             , 0.0247;
       'csv_files\30G.csv'             , 0.0297;
       'csv_files\35G.csv'             , 0.0345;
       'csv_files\50.2 G.csv'          , 0.0502;
       'csv_files\55G.csv'             , 0.055;
       'csv_files\60G.csv'             , 0.060;
       'csv_files\64.8G.csv'           , 0.0648;
       'csv_files\70.1G.csv'           , 0.0701;
       'csv_files\84 G part 2.csv'    , 0.0847;   
    };

%% code:
% Fit each result:

CycleTime = zeros(1,size(lab_results,1));
damp_style = 'clean_';
if add_linear_damp
    damp_style = 'lin_';
end 

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
        
        %Fix zero
        graph_zero = y(length(y));
        if abs(graph_zero) > distance_error_range
            for j = 1:length(y)
                y(j) = y(j) - graph_zero;
            end
        end
        
        
        %% Remove unwanted data
        
        % Remove data captured before the experiment begins:
        if find(y == max(y)) < find(y == min(y))
            init_wave = find(y == max(y));
        else
            init_wave = find(y == max(y));
        end
        x = x-x(init_wave);
        for j = 1:init_wave-3
            x(1) = [];
            y(1) = [];
        end
        
        %Remove data captured after the experiment ends:
        for j = 60:length(x)
            if (abs(y(j))<distance_error_range) && (abs(y(j-15))<distance_error_range)
                x(j:length(x)) = [];
                y(j:length(y)) = [];
            break
            end
        end
        
        %Remove data repeatitions
        [b,m1,n1] = unique(x,'first');
        [c1,d1] =sort(m1);
        x = b(d1);
        y = y(m1);
        final_fit = DampedHarmonic_fit(x,y, add_linear_damp);
        
        % Find peaks:
        [pks,locs] = findpeaks(y,x);
        amplitude_fit = DampedAmplitudeFit(locs,pks, add_linear_damp);
        peak_error = zeros(1,length(locs)) + distance_error_range;
        time_error = zeros(1,length(locs));
        %CycleTime(i) = mean(diff(locs(2:10)));
        fit_values = coeffvalues(final_fit);
        CycleTime(i) = fit_values(3);
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
            f = gcf;
            exportgraphics(f,[image_save_path 'part_1_' damp_style char(string(lab_results(i,2))) '.png'],'Resolution',300);

            %Plot peaks per time
            figure
            hold on
            errorbar(locs, pks, peak_error, peak_error, time_error, time_error, 'color','magenta','LineStyle','none', 'LineWidth', 2)
            plot(amplitude_fit)
            grid
            box on
            xlabel('Time(S)')
            ylabel('peaks(M)')
            legend('peaks', 'Fitted Curve')
            if ~add_linear_damp
                set(gca,'YScale','log');
            end
            hold off
            f = gcf;
            exportgraphics(f,[image_save_path 'part_2_' damp_style char(string(lab_results(i,2))) '.png'],'Resolution',300);
            
        end
    end
end

%% Plot for all mass
%figure
mass = cell2mat(lab_results(:,2))'.*1000;
omega = 2*pi./CycleTime;
%plot(mass,omega, '.')
hold on
xlabel('Mass(KG)')
ylabel('omega(Rad/S)')
legend('peaks', 'Fitted Curve')
hold off

function f = DampedHarmonic_fit(x, y, add_linear_damp)
    %% Get fit parametes:
    
    y = detrend(y);                                                                      % Remove Linear Trend
    yu = max(y);
    yl = min(y);
    yr = (yu-yl);                                                                        % Range of  y’
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                                 % Returns Approximate Zero-Crossing Indices Of Argument Vector
    zt = x(zci(y));
    per = 2*mean(diff(zt));                                                              % Estimate period
    ym = mean(y);                                                                        % Estimate offset
    
    if add_linear_damp
        init_fit = @(b,x)  b(1) .* exp(b(2).*x-b(3)*x) .* (cos(2*pi*x./b(4) + 2*pi/b(5))) + b(6);   % Objective Function to fit
        fcn = @(b) norm(init_fit(b,x) - y);                                              % Least-Squares cost function
        [s,] = fminsearch(fcn, [yr; -10; 0;  per;  -1;  ym]);                                   % Minimise Least-Squares
        fit_params = s';
        %% Fit:
        fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', fit_params);         % Use the parameters gathered as starting points.
        fitt = fittype('(a.*exp(b*x)- c*x).*(cos(2*pi.*x/d + 2.*pi./e)) + f','coefficients', {'a', 'b', 'c', 'd', 'e', 'f'}, 'options', fo);
        f = fit(x,y,fitt);    
    else
        init_fit = @(b,x)  b(1) .* exp(b(2).*x) .* (cos(2*pi*x./b(3) + 2*pi/b(4))) + b(5);   % Objective Function to fit
        fcn = @(b) norm(init_fit(b,x) - y);                                              % Least-Squares cost function
        [s,] = fminsearch(fcn, [yr; -10; per;  -1;  ym]);                                   % Minimise Least-Squares        
        fit_params = s';
        %% Fit:
        fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', fit_params);         % Use the parameters gathered as starting points.
        fitt = fittype('a.*exp(b*x).*(cos(2*pi.*x/c + 2.*pi./d)) + e','coefficients', {'a', 'b', 'c', 'd', 'e'}, 'options', fo);
        f = fit(x,y,fitt);
    end
end

function f = DampedAmplitudeFit(x, y, add_linear_damp)
    linear_damp = '+ c*x';
    coefficients = {'a', 'b', 'c'};
    if ~add_linear_damp
        linear_damp = '';
        coefficients = {'a', 'b'};
    end
    fitt = fittype(['a*exp(-b*x)' linear_damp],'coefficients', coefficients);
    f = fit(x,y,fitt);
end

