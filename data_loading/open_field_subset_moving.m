function in_motion = open_field_subset_moving(ds)
if ds.num_trials ~= 1
    error('must be open field data with 1 "trial"');
end
XY = ds.trials.centroids;
range = max(XY) - min(XY);
cmppx = 1/(range/[46 36]); %should include box dimensions in cm

idt = 20;
V = cmppx*diff(XY)*idt; %in cm/s
avg_over = 1; %one second
Vd = medfilt1(V, avg_over*idt);
speed_d = hypot(Vd(:,1),Vd(:,2));

tHi = 1;
tLo = 0.3;

in_motion = h_th(speed_d, tHi, tLo);
end

function mov = h_th(s, tHi, tLo)
mov = false(length(s) + 1,1);
moving = false;
for i = 1:length(s)
    if moving
        moving = s(i) > tLo;
    else
        moving = s(i) > tHi;
    end
    mov(i) = moving;
end
mov(end) = mov(end-1);
end