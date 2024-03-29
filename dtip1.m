function output_txt = dtip_function(obj,event_obj,tBld,tSkl,oBld,oSkl)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
dindex = get(event_obj,'DataIndex');

output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)],['I:',num2str(dindex)],...
    ['T-BLD:',num2str(tBld(dindex))]};


% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end
