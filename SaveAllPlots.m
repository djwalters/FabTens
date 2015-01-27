hfigs = get(0,'children');  %Get list of figures

LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
directory = uigetdir(LocalPath,...
            'Select  directory to save plots');

for i = 1:length(hfigs)    
    figure(hfigs(i))        %Bring figure to foreground
    filename = get(gcf,'name'); %Get window title for filename
    FilePath = fullfile(directory,filename);
    saveas(hfigs(i),FilePath,'fig');
    saveas(hfigs(i),FilePath,'emf');
    saveas(hfigs(i),FilePath,'png');
end
