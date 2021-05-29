function encoder_data = convert_encoder(data,timestamps)
dt = timestamps(2)-timestamps(1);
try
rotate_pulse = data';
catch ME
    disp('No ADC detected')
    return
end
pwd_max = max(rotate_pulse);
pwd_min = min(rotate_pulse);
pwd_threshold = pwd_max/2;
i = 1;
T(:,1) = (1:length(rotate_pulse))*dt;
pulse_detection{i} = [];
    for ii = 2:length(rotate_pulse)
        if abs(rotate_pulse(ii,i))>pwd_threshold && abs(rotate_pulse(ii-1,i))<=pwd_threshold
            pulse_detection{i} = [pulse_detection{i}, ii*dt]; 
        end
    end
pulse_detectionf(:,2) = (1:length(pulse_detection{1}));
pulse_detectionf(:,1) = cell2mat(pulse_detection)';
% dfpulse_detection(:,2) = (diff(pulse_detectionf(:,2)));
% dfpulse_detection(:,1) = (1:length(dfpulse_detection));
%1000 pulse detections per linar rotation gives us 0.360 degrees of change
radius = 7.62; %radius of disk (cm)
theta = 0.360; theta_total = theta.*pulse_detectionf(:,2);
time = pulse_detectionf(:,1);
ang_distance = theta_total.*radius;
velocity_theta = diff(theta_total);
velocity_time = diff(time);
radian_avg = velocity_theta./velocity_time;
ang_velocity = radian_avg.*radius;

encoder_data.dt = dt;
encoder_data.rotate_pulse = rotate_pulse;
encoder_data.pulse_detectionf = pulse_detectionf;
encoder_data.ang_distance = ang_distance;
encoder_data.velocity_theta = velocity_theta;
encoder_data.velocity_time = velocity_time;
encoder_data.ang_velocity = ang_velocity;
encoder_data.time = time;
end
