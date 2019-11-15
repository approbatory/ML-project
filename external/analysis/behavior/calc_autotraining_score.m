function [ percent_correct ] = calc_autotraining_score( file, correct_turn_direction )
%calc_autotraining_score Reads autotraining text file and gives behavior
%file should be a string
%correct_turn_direction should be a string, e.g. 'left' or 'right'
%score as % correct
%2014-12-17 Fori Wang

    fid = fopen(file);
    maze_data = textscan(fid,'%s %s %s %f');
    fclose(fid);
    
    actual_turn_direction = maze_data{3};
    num_trials = length(actual_turn_direction);
    
    num_correct_instances = sum(strcmp(correct_turn_direction,actual_turn_direction));
    percent_correct = num_correct_instances*100/num_trials;

end