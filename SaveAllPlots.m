% Script saves all open figure windows to 3 file formats to the directory
% specified by the user.  (3 formats: .fig, .emf, .png)
hfigs = get(0,'children');  %Get list of figures

LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
directory = uigetdir(LocalPath,...
            'Select  directory to save plots');

for i = 1:length(hfigs)
    % Set figure size and position on screen
    set(hfigs(i),'Units','inches')
    set(hfigs(i),'Position',[3,3,6,3.5])
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
    saveas(hfigs(i),FilePath,'emf');
    saveas(hfigs(i),FilePath,'png');
end
