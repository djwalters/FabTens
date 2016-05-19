% Script saves all open figure windows to 3 file formats to the directory
% specified by the user.  (3 formats: .fig, .emf, .png)
hfigs = get(0,'children');  %Get list of figures

LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
directory = uigetdir(LocalPath,...
    'Select  directory to save plots');
DataPath = fullfile(directory,'TensorData');
save(DataPath)

for i = 1:length(hfigs)
    % Set figure size and position on screen
    set(hfigs(i),'PaperUnits','inches')
    set(hfigs(i),'PaperPosition',[3,3,6,3.5])
    set(hfigs(i),'PaperPositionMode','manual')
end
for i = 1:length(hfigs)
    figure(hfigs(i))        %Bring figure to foreground
    filename = get(gcf,'name'); %Get window title for filename
    % If no special window title specified, give generic filename
    if isempty(filename)
        fname = sprintf('Figure %2.0f',i);
        filename = fname;
    end
    FilePath = fullfile(directory,filename);
    % Save figures to specified path
    saveas(hfigs(i),FilePath,'fig');
    print(hfigs(i),FilePath,'-dmeta');
    print(hfigs(i),FilePath,'-dpng');
end
