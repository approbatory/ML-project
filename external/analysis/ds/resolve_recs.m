function res_list = resolve_recs(md)

h = figure;

% For matched cell index k, "resolution_list" indicates:
%   res_list(k,1): Selected rec index
%   res_list(k,2): Cell index within that rec
res_list = zeros(md.num_cells, 2);

% If filter only shows up on one rec, then assign it automatically
for j = 1:md.num_cells
    inds = md.matched_indices(j,:);
    num_recs = sum(inds>0);
    if (num_recs == 1)
        rec_ind = find(inds, 1);
        res_list(j,:) = [rec_ind inds(rec_ind)];
    end
end

cell_idx = 1;
while (1)
    draw_cell(cell_idx);
    
    % Ask user for command
    if (res_list(cell_idx) ~= 0)
        assign_str = sprintf('Rec %d', res_list(cell_idx,1));
    else
        assign_str = 'unassigned';
    end
    prompt = sprintf('Resolve recs (ID %d of %d: %s) >> ',...
        cell_idx, md.num_cells, assign_str);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number
        if ((1 <= val) && (val <= md.num_days)) % Check valid rec index
            if (md.matched_indices(cell_idx, val) ~= 0)
                res_list(cell_idx,1) = val;
                res_list(cell_idx,2) = ...
                    md.matched_indices(cell_idx, val);

                go_to_next_cell();
            else
                fprintf('  Sorry, %d is not a valid rec index for this cell\n', val);
            end
        else
            fprintf('  Sorry, %d is not a valid rec index\n', val);
        end
    else
        resp = lower(resp);
        if isempty(resp) % Empty string gets mapped to "n"
            resp = 'n';
        end
        
        switch resp(1)
            case 'f'
                cell_idx = 1;

            case 'c' % Jump to a cell
                val = str2double(resp(2:end));
                if ~isnan(val)
                    if ((1 <= val) && (val <= md.num_cells))
                        cell_idx = val;
                    else
                        fprintf('  Sorry, %d is not a valid cell on Day %d\n', val, md.sort_day);
                    end
                end
                
            case 'n'
                go_to_next_cell();
                
            case 'p' % Previous cell
                if (cell_idx > 1)
                    cell_idx = cell_idx - 1;
                else
                    fprintf('  Already at first md cell!\n');
                end
                
            case 'q' % Exit
                close(h);
                break;
                
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end    
end

    function go_to_next_cell()
        to_be_assigned = (res_list(:,1) == 0);
        next_idx = find_next_cell_to_process(cell_idx, to_be_assigned);
        if isempty(next_idx)
            fprintf('  All cells have been assigned!\n');
        else
            cell_idx = next_idx;
        end
    end

    function draw_cell(common_cell_idx)
        clf;
        
        % Draw filters
        for k = 1:md.num_days
            day = md.valid_days(k);
            cell_idx_k = md.get_cell_idx(common_cell_idx, day);
            
            if (cell_idx_k ~= 0)
                % Draw filters
                ds_cell = md.day(day).cells(cell_idx_k);
                subplot(3, md.num_days, k);
                com = ds_cell.com;
                boundary = ds_cell.boundary;
                
                imagesc(ds_cell.im);
                colormap gray;
                axis image;
                hold on;
                plot(boundary(:,1),boundary(:,2),'Color',get_color(k),...
                     'LineWidth', 2);
                hold off;
                zoom_half_width = min(size(ds_cell.im))/20;
                xlim(com(1)+zoom_half_width*[-1 1]);
                ylim(com(2)+zoom_half_width*[-1 1]);
                
                title_str = sprintf('Iter %d -- Cell %d', day, cell_idx_k);
                title(title_str);
            end
        end
        
        % Draw traces       
        h_trace = subplot(3, md.num_days, (md.num_days+1):(3*md.num_days));
        set(zoom(h_trace), 'Motion', 'horizontal');
        trace_offset = 0;
        for k = fliplr(1:md.num_days)
            day = md.valid_days(k);
            cell_idx_k = md.get_cell_idx(common_cell_idx, day);
            
            if (cell_idx_k ~= 0)
                trace_k = md.day(day).get_trace(cell_idx_k);
                plot(trace_k + trace_offset, 'Color', get_color(k));
                trace_offset = trace_offset + max(trace_k);
                hold on;
            end
        end
        xlim([0 length(trace_k)]);
        ylim([0 trace_offset]);
        grid on;
        xlabel('Frame');
        ylabel('Traces');
        set(gca, 'YTickLabel', []);

        function color = get_color(day_idx)
            colors = {'b', 'r', [0 0.5 0]};
            color = colors{mod(day_idx-1, length(colors))+1};
        end 
        
    end % draw_cell

end % resolve_recs