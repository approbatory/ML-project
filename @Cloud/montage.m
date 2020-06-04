function montage
figure;
idx = 1;
s_ = @(i)subplot(3,4,i);
Mouse_origins = {'Mouse2022', 'Mouse2024', 'Mouse2028'};
C{1} = Cloud(DecodeTensor.special_by_mouse(Mouse_origins{1}));
C{2} = Cloud(DecodeTensor.special_by_mouse(Mouse_origins{2}));
C{3} = Cloud(DecodeTensor.special_by_mouse(Mouse_origins{3}));
for i = 1:3
    s_(idx);
    C{i}.cum_loadings;
    idx = idx + 1;
    s_(idx);
    C{i}.signal_geo;
    idx = idx + 1;
    s_(idx);
    C{i}.noise_geo;
    idx = idx + 1;
    s_(idx);
    C{i}.signal_noise_overlap_geo;
    idx = idx + 1;
end
disp(Mouse_origins);
end