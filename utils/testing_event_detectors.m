load ../../linear_track/Mouse2022/Mouse-2022-20150326-linear-track/Mouse-2022-20150326_093722-linear-track-TracesAndEvents.mat

%%

%choosing cell 345
c_i = 345;

fprintf('getting original trace:'); tic
dff_trace = tracesEvents.rawTraces(:,c_i);
fprintf(' done\n'); toc

fprintf('getting spikeDeconv version:'); tic
pj_events = tracesEvents.spikeDeconv(:,c_i);
fprintf(' done\n'); toc

%fprintf('getting mlspike calibrated version:'); tic
%freq = 20;
%mls_event_times = mlspike_yescal(dff_trace, 1/freq);
%mls_events = 0*dff_trace;
%mls_events(round(mls_event_times*freq)) = 1;
%fprintf(' done\n'); toc

fprintf('getting computeSignalPeaks version:'); tic
csp_events = computeSignalPeaks(dff_trace')';
fprintf(' done\n'); toc
%%
fprintf('getting onetimeevents version:'); tic
ote_events = one_time_events(dff_trace')';
fprintf(' done\n'); toc

fprintf('getting FST version:'); tic
fst_events = Utils.event_detection(dff_trace, false);
fprintf(' done\n'); toc

fprintf('getting IED version:'); tic
ied_events = iterative_event_detection(dff_trace);
fprintf(' done\n'); toc
%%

figure;
hold on;

num_frames = numel(dff_trace);
time = (1:num_frames)/freq;

plot(time,dff_trace, 'DisplayName', 'DFF');
plot(time,pj_events~=0, '-o', 'DisplayName', 'spikeDeconv');
%plot(time,mls_events, 'DisplayName', 'MLspike');
%plot(time,csp_events, '-o', 'DisplayName', 'computeSignalPeaks');
%plot(time,-2*ote_events, 'DisplayName', 'OTE');
%plot(time, 2*fst_events, 'DisplayName', 'FST');
plot(time, 2-4*ied_events, 'DisplayName', 'IED');
tracesEvents.position(1:91,1)=600;
plot(time, -2+tracesEvents.position(:,1)./range(tracesEvents.position(:,1)), '.', 'DisplayName', 'position');

legend;