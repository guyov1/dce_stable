function mouseMoveOverDataWindow(hObject,~)

handles=guidata(hObject);
pos = get(hObject,'CurrentPoint');
% disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);
% handles.mouse_pos=pos;
hObject.mouse_pos=pos;
guidata(hObject,handles);
end