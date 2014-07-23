function output_txt = dtip_function(obj,event_obj,tBld,tSkl,oBld,oSkl)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
%dindex = get(event_obj,'DataIndex');
dindex = 1;
output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)],['I:',num2str(dindex)],...
    ['tBLD',num2str(tBld(dindex))],...
    ['tSKL',num2str(tSkl(dindex))],...
    ['jBLD',num2str(oBld(dindex))],...
    ['jSKL',num2str(oSkl(dindex))]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end
