function c = get_cluster(partition, memory, time, cpus, extra)

if ~exist('extra', 'var')
    extra = '';
end

addpath('/oak/stanford/groups/mschnitz/bahanonu_shared/mathworks/matlab')

configCluster;

% Load the cluster and setup
c = parcluster;
% Int: number of processes per node
c.AdditionalProperties.ProcsPerNode = cpus;
% Str: which partition to run on
c.AdditionalProperties.QueueName = partition;
% Submit to normal and H&S partitions, DO NOT use mschnitz
c.AdditionalProperties.AdditionalSubmitArgs = extra;
% How long the pool can be open
c.AdditionalProperties.WallTime = time;
% What address to email upon completion or error
c.AdditionalProperties.EmailAddress = 'omer.hazon@gmail.com';
% How much memory to request per worker
c.AdditionalProperties.MemUsage = memory;
c.AdditionalProperties.DebugMessagesTurnedOn = true;
c.NumThreads = 1;
c.saveProfile;
%c.AdditionalProperties