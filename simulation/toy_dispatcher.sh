S=8
echo "calling $S matlabs"
for i in $(seq $S); do
    matlab -nodesktop -nojvm -r 'Q = pwd; cd ~; addpath(pathdef); cd(Q); toy_sampler; exit' &
done

echo "waiting for them to finish"
wait

echo "done"
