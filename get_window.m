function window = get_window(data, idx, width)
    % Return mean of data over a window of specified width around a
    % specified point
    window = mean(data(idx-width:idx+width));
end

