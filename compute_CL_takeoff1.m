function compute_CL_takeoff()
    % INPUT PARAMETERS
    mass = 3;                         % kg
    thrust = 30;                      % N (from APC datasheet)
    rho = 1.225;                      % kg/m^3 (air density at sea level)
    s = input('Enter lifting surface area in m^2: ');  % e.g. 0.5 m^2
    aoa_deg = input('Enter AoA range as [min max] in degrees: '); % e.g. [10 30]
    
    % ANGLE OF ATTACK RANGE (in radians for calculations)
    aoa_range = linspace(aoa_deg(1), aoa_deg(2), 100);
    aoa_rad = deg2rad(aoa_range);
    
    % CONSTANTS
    g = 9.81;     % m/s^2
    W = mass * g; % Weight (N)
    takeoff_distance = 3;  % meters

    % INITIALIZE OUTPUT
    required_CL = zeros(size(aoa_range));
    takeoff_speed = zeros(size(aoa_range));
    
    for i = 1:length(aoa_range)
        alpha = aoa_rad(i);
        
        % Thrust components
        T_x = thrust * cos(alpha);  % Forward (acceleration)
        T_y = thrust * sin(alpha);  % Vertical (lift assist)

        % Net required lift from wings
        L_required = W - T_y;

        % Acceleration (a = F/m)
        accel = T_x / mass;

        % Using kinematic equation: s = (v^2)/(2a) => solve for v
        V = sqrt(2 * accel * takeoff_distance);
        takeoff_speed(i) = V;

        % Now use lift equation to find CL
        CL = (2 * L_required) / (rho * V^2 * s);
        required_CL(i) = CL;
    end

    % PLOT RESULTS
    figure;
    plot(aoa_range, required_CL, 'LineWidth', 2);
    grid on;
    xlabel('Angle of Attack (degrees)');
    ylabel('Required C_L');
    title('Required Coefficient of Lift vs. Angle of Attack');
    
    % OPTIONAL: Display table of values
    fprintf('\n%-10s %-15s %-15s\n','AoA (deg)','CL Required','Takeoff Speed (m/s)');
    for i = 1:5:length(aoa_range)
        fprintf('%-10.2f %-15.4f %-15.4f\n', aoa_range(i), required_CL(i), takeoff_speed(i));
    end
end