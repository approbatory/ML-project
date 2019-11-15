function [ percent_correct ] = calc_behavior_score( file )
% Reads behavior text file and gives behavior score as % correct
%
% Updated 2015-03-12 Fori Wang

    [~,maze_data,~] = parse_plusmaze(file);
    end_arm = maze_data(:,2);
    actual_end_arm = maze_data(:,3);
    num_trials = size(maze_data,1);
    
    num_correct_instances = sum(strcmp(end_arm,actual_end_arm));
    percent_correct = num_correct_instances*100/num_trials;

end