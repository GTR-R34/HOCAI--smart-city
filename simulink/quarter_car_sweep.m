% Quarter-car sweep -> labeled training CSV
rng(1);
fs_sim = 1000; fs_out = 100;                 % simulate at 1 kHz, sensor at 100 Hz
t = (0:1/fs_sim:2-1/fs_sim)';                % 2 s window
g = 9.81; dt = 1/fs_sim;
classes = {'smooth','rough','speed_bump','pothole'};
N = 300;
out = cell(numel(classes)*N,1); runId = 0;

for c = 1:numel(classes)
    for k = 1:N
        runId = runId + 1;
        % vehicle parameters (+-20% mass)
        ms = 250*(0.8+0.4*rand); mu = 40; ks = 18000; cs = 1000; kt = 200000;
        A = [0 1 0 0; -ks/ms -cs/ms ks/ms cs/ms; 0 0 0 1; ks/mu cs/mu -(ks+kt)/mu -cs/mu];
        B = [0;0;0;kt/mu];  C = [-ks/ms -cs/ms ks/ms cs/ms];
        Ad = expm(A*dt);  Bd = A\(Ad-eye(4))*B;

        spd = 10+50*rand;  v = spd/3.6;          % km/h and m/s
        z = zeros(size(t));  t0 = 0.6+0.4*rand;
        depth = 0; width = 0;
        switch classes{c}
            case 'pothole'
                depth = 0.02+0.13*rand;  width = 0.2+0.6*rand;
                d = width/v;  m = t>=t0 & t<t0+d;
                z(m) = -depth*sin(pi*(t(m)-t0)/d);
            case 'speed_bump'
                depth = 0.04+0.06*rand;  width = 1+2*rand;
                d = width/v;  m = t>=t0 & t<t0+d;
                z(m) = depth*sin(pi*(t(m)-t0)/d);
            case 'rough'
                z = movmean(randn(size(t)),20)*(0.01+0.02*rand);
            case 'smooth'
                z = movmean(randn(size(t)),20)*0.001;
        end

        x = zeros(4,1); acc = zeros(size(t));
        for i = 1:numel(t)
            acc(i) = C*x;
            x = Ad*x + Bd*z(i);
        end
        acc = movmean(acc,10);                    % simple anti-alias
        idx = 1:fs_sim/fs_out:numel(t);
        a = acc(idx)/g;  n = numel(idx);
        az = 1 + a + 0.02*randn(n,1);
        ax = 0.1*a + 0.02*randn(n,1);
        ay = 0.02*randn(n,1);
        out{runId} = table(repmat(runId,n,1), (0:n-1)'/fs_out, ax, ay, az, ...
            repmat(string(classes{c}),n,1), repmat(round(spd),n,1), ...
            repmat(round(depth*100,1),n,1), repmat(round(width*100),n,1), ...
            'VariableNames',{'run_id','t','ax','ay','az','label','speed_kmh','depth_cm','width_cm'});
    end
end
writetable(vertcat(out{:}),'training_data.csv');
disp('done: training_data.csv');