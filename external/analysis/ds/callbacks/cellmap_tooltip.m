function tooltip_txt = cellmap_tooltip(~, event_obj)
    pos = get(event_obj, 'Position');
    tooltip_txt = sprintf('ID: %d', pos(3));
end