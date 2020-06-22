function output_txt = update_data_tip_fn(Settings, obj, event_obj)
% Display the position of the data cursor
%
% Arguments:
% -----------------------
% Settings: settings class object
%
% obj          Currently not used (empty)
%
% event_obj    Handle to event object
%
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj, 'Position');
index = get(event_obj, 'DataIndex');
hAxes = gca;
mCurrent = getimage(hAxes);
value = mCurrent(index);

% get global coordinates
ROI = Settings.WindowLimits;
x_local = pos(1);
y_local = pos(2);
xform_params = [ROI(1), ROI(3), Settings.SoftwareBinSizeH, Settings.SoftwareBinSizeV];
[x_global, y_lgobal] = roi2img_coord(xform_params, x_local, y_local);

output_txt = {['X full img: ',num2str(x_global, 8)],...
              ['Y full img: ',num2str(y_lgobal, 8)],...
              ['X roi: ', num2str(x_local, 8)],...
              ['Y roi: ', num2str(y_local, 8)],...
              ['Index: ',num2str(value, 6)]};
end