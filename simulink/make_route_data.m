% make_route_data.m
% Creates: route_telemetry.csv (live feed), route_truth.csv (answer key),
%          route_segments.csv (road coordinates)
rng(2);

%% 1. SETTINGS
fs_sim = 1000; fs_out = 100; g = 9.81; dt = 1/fs_sim;
base_ts = 1791800000;      % start time in seconds (any number works)
nDays = 10; nVeh = 5;

% Road segments: [lat1 lon1 lat2 lon2]
seg = [12.9352 77.5345 12.9370 77.5345;    % Road 1 (ambulance's best road, damages fast)
       12.9352 77.5365 12.9370 77.5365;    % Road 2 (alternative, damages slowly)
       12.9340 77.5345 12.9340 77.5363];   % Road 3 (healthy, has a speed bump)
nSeg = size(seg,1);
seg_len = zeros(nSeg,1);                   % length in metres
for s = 1:nSeg
    dlat = (seg(s,3)-seg(s,1))*111320;
    dlon = (seg(s,4)-seg(s,2))*111320*cosd(seg(s,1));
    seg_len(s) = sqrt(dlat^2 + dlon^2);
end

% Vehicles
veh_id    = ["bus_01","bus_02","bus_03","car_04","van_05"];
veh_speed = [25 30 35 40 30];              % km/h
veh_mass  = 250*(0.8+0.4*rand(1,nVeh));    % each vehicle has its own body mass

% Hazards: [segment type pos_m appear_day size_cm width_m]   type: 1=pothole, 2=speed bump
H = [ 1 1  30 2  6 0.4;
      1 1  70 3  8 0.5;
      1 1 100 4  5 0.4;
      1 1 130 6  9 0.5;
      1 1 160 7  7 0.4;
      1 1 185 9 10 0.6;
      2 1  60 5  4 0.4;
      2 1 140 9  5 0.4;
      3 2 100 1  7 2.0 ];

%% 2. DRIVE EVERY PASS
out = cell(nDays*nVeh*nSeg,1); n = 0;
for day = 1:nDays
  for v = 1:nVeh
    for s = 1:nSeg
      n = n + 1;
      L   = seg_len(s);
      spd = max(veh_speed(v) + 3*randn, 15);   % km/h with small variation
      vms = spd/3.6;                           % m/s
      T   = L/vms;                             % time to drive the segment
      t   = (0:dt:T)';

      % road height: tiny roughness + every hazard that exists by today
      z = movmean(randn(size(t)),20)*0.002;
      for h = find(H(:,1)==s & H(:,4)<=day)'
          typ = H(h,2);  size_m = H(h,5)/100;
          if typ == 1
              size_m = size_m*(1 + 0.05*(day-H(h,4)));   % potholes grow 5% per day
              sgn = -1;                                  % dip
          else
              sgn = +1;                                  % hump
          end
          t0 = H(h,3)/vms;  d = H(h,6)/vms;
          mk = t>=t0 & t<t0+d;
          z(mk) = z(mk) + sgn*size_m*sin(pi*(t(mk)-t0)/d);
      end

      % physics -> sensor readings
      acc = simulate_pass(veh_mass(v), z, dt);
      acc = movmean(acc,10);
      idx = 1:fs_sim/fs_out:numel(t);
      a   = acc(idx)/g;   k = numel(idx);
      az  = 1 + a + 0.02*randn(k,1);
      ax  = 0.1*a + 0.02*randn(k,1);
      ay  = 0.02*randn(k,1);

      % GPS and timestamp
      tt   = t(idx);
      frac = min(vms*tt/L, 1);
      lat  = seg(s,1) + frac*(seg(s,3)-seg(s,1)) + 1e-5*randn(k,1);
      lon  = seg(s,2) + frac*(seg(s,4)-seg(s,2)) + 1e-5*randn(k,1);
      start_off = 8*3600 + (v-1)*600 + (s-1)*120 + 30*rand;
      ts = base_ts + (day-1)*86400 + start_off + tt;

      out{n} = table(repmat(veh_id(v),k,1), ts, lat, lon, ax, ay, az, ...
          'VariableNames',{'vehicle_id','ts','lat','lon','ax','ay','az'});
    end
  end
end

%% 3. SAVE THE LIVE FEED
T = sortrows(vertcat(out{:}),'ts');
writetable(T,'route_telemetry.csv');

%% 4. SAVE THE ANSWER KEY
nH   = size(H,1);
type = repmat("pothole",nH,1);  type(H(:,2)==2) = "speed_bump";
frac = H(:,3)./seg_len(H(:,1));
tlat = seg(H(:,1),1) + frac.*(seg(H(:,1),3)-seg(H(:,1),1));
tlon = seg(H(:,1),2) + frac.*(seg(H(:,1),4)-seg(H(:,1),2));
truth = table((1:nH)', H(:,1), type, H(:,3), tlat, tlon, H(:,4), H(:,5), ...
    'VariableNames',{'hazard_id','segment','type','pos_m','lat','lon','appear_day','size_cm'});
writetable(truth,'route_truth.csv');

%% 5. SAVE THE ROAD COORDINATES
segtab = table((1:nSeg)', seg(:,1), seg(:,2), seg(:,3), seg(:,4), seg_len, ...
    'VariableNames',{'segment','lat1','lon1','lat2','lon2','length_m'});
writetable(segtab,'route_segments.csv');

fprintf('done: %d rows in route_telemetry.csv\n', height(T));

%% ---- local function (must stay at the bottom of the file) ----
function acc = simulate_pass(ms, z, dt)
    mu = 40; ks = 18000; cs = 1000; kt = 200000;
    A = [0 1 0 0; -ks/ms -cs/ms ks/ms cs/ms; 0 0 0 1; ks/mu cs/mu -(ks+kt)/mu -cs/mu];
    B = [0;0;0;kt/mu];  C = [-ks/ms -cs/ms ks/ms cs/ms];
    Ad = expm(A*dt);  Bd = A\(Ad-eye(4))*B;
    x = zeros(4,1);  acc = zeros(size(z));
    for i = 1:numel(z)
        acc(i) = C*x;
        x = Ad*x + Bd*z(i);
    end
end