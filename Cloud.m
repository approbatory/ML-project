classdef Cloud

    properties
        dt;
        clouds;
        clouds_shuf;
        mus;
        evecs;
        evecs_shuf;
        dmus;
        loadings;
        loadings_shuf;
        K;
        N = 21;
    end
    
    methods
        function o = Cloud(i)
            if isa(i, 'DecodeTensor')
                o.dt = i;
            else
                o.dt = DecodeTensor.cons_filt(i);
            end
            
            n = @(x)x./norm(x);
            
            [N, o.K, T] = size(o.dt.data_tensor);
            for j = 1:2*o.K
                o.clouds{j} = o.dt.get_bin_resp(j);
                o.clouds_shuf{j} = Cloud.shuffle(o.clouds{j});
                o.mus{j} = mean(o.clouds{j},2);
                o.evecs{j} = pca(o.clouds{j}');
                o.evecs_shuf{j} = pca(o.clouds_shuf{j}');
            end
            
            for j = 1:o.K-1
                o.dmus{j} = o.mus{j+1} - o.mus{j};
                o.loadings{j} = abs(n(o.dmus{j})' * o.evecs{j});
                o.loadings_shuf{j} = abs(n(o.dmus{j})' * o.evecs_shuf{j});
            end
            for j = o.K+1:2*o.K-1
                o.dmus{j} = o.mus{j+1} - o.mus{j};
                o.loadings{j} = abs(n(o.dmus{j})' * o.evecs{j});
                o.loadings_shuf{j} = abs(n(o.dmus{j})' * o.evecs_shuf{j});
            end
            
            o.N = o.cum_loadings(true);
        end
        
        function N_50 = cum_loadings(o, suppress) %1
            if ~exist('suppress', 'var')
                suppress = false;
            end
            
            hold on
            %colors = hsv(o.K);
            for i = 1:2*o.K-1
                if isempty(o.loadings{i})
                    continue;
                end
                cs = cumsum(o.loadings{i}.^2);
                cs_s = cumsum(o.loadings_shuf{i}.^2);
                reorder = @(x) x(randperm(numel(x)));
                cs_ro = cumsum(reorder(o.loadings{i}.^2));
                %plot(cs, 'Color', colors(mod(i-1,o.K)+1,:));
                if ~suppress
                    plot(cs, 'b');
                    plot(cs_s, 'r');
                    plot(cs_ro, 'g');
                end
                ind(i) = find(cs > 0.5, 1);
            end
            ind = ind(ind ~= 0);
            %figure;
            %histogram(ind);
            N_50 = round(median(ind));
            fprintf('Median index to 50%% is %d\n', N_50);
            
            if ~suppress
                xlabel 'Number of PCs'
                ylabel 'Fraction of \Delta\mu within the subspace'
                text(160, 0.2, 'Real', 'Color', 'b');
                text(160, 0.15, 'Shuffled', 'Color', 'r');
                text(160, 0.1, 'Unordered', 'Color', 'g');
                title 'Showing all spatial bins'
            end
        end
        
        function signal_geo(o) %2
            n = @(x)x./norm(x);
            dm = cellfun(n, o.dmus, 'UniformOutput', false);
            
            cos_overlap = zeros(2*o.K-1);
            for i = 1:2*o.K-1
                if isempty(dm{i}), continue; end
                for j = 1:2*o.K-1
                    if isempty(dm{j}), continue; end
                    cos_overlap(i,j) = dm{i}'*dm{j};
                end
            end
            
            
            imagesc(cos_overlap);
            xlabel 'Spatial bin'
            ylabel 'Spatial bin'
            colormap(gca, bluewhitered);
            colorbar;
            title 'Cos overlap between signal directions'
            axis image
        end
        
        function noise_geo(o) %3
            cos_overlap = zeros(2*o.K);
            for i = 1:2*o.K
                for j = 1:2*o.K
                    cos_overlap(i,j) = mean(cos(subspacea(o.evecs{i}(:,1:o.N), o.evecs{j}(:,1:o.N))));
                end
            end
            
            
            imagesc(cos_overlap);
            line([20 20]+0.5, ylim, 'Color', 'w');
            line(xlim, [20 20]+0.5, 'Color', 'w');
            xlabel 'Spatial bin';
            ylabel 'Spatial bin';
            colorbar;
            title(sprintf('Mean canon corr between\n%dD noise subspaces', o.N));
            axis image
        end
        
        function signal_noise_overlap_geo(o) %4
            n = @(x)x./norm(x);
            dm = cellfun(n, o.dmus, 'UniformOutput', false);
            
            
            cos_overlap = zeros(2*o.K, 2*o.K-1);
            for i = 1:2*o.K
                for j = 1:2*o.K-1
                    if isempty(dm{j})
                        cos_overlap(i,j) = nan;
                        continue;
                    end
                    cos_overlap(i,j) = cos(subspace(o.evecs{i}(:,1:o.N), dm{j}));
                end
            end
            
            
            imagesc(cos_overlap);
            line([20 20]+0.5, ylim, 'Color', 'w');
            line(xlim, [20 20]+0.5, 'Color', 'w');
            xlabel 'Spatial bin of signal direction';
            ylabel 'Spatial bin of noise subspace';
            colorbar;
            title(sprintf('Cos overlap between %dD noise subspaces\nand signal directions', o.N))
            axis image;
        end
        
        function summary(o)
            figure;
            subplot(1,4,1);
            o.cum_loadings;
            subplot(1,4,2);
            o.signal_geo;
            subplot(1,4,3);
            o.noise_geo;
            subplot(1,4,4);
            o.signal_noise_overlap_geo;
        end
    end
    
    methods(Static)
        function X = shuffle(X)
            [n, m] = size(X);
            for i = 1:n
                X(i,:) = X(i, randperm(m));
            end
        end
        
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
    end

end