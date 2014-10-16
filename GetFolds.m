function TimeFolds = GetFolds(RootPath)

d = dir(RootPath);
isub = [d(:).isdir]; %# returns logical vector
TimeFolds = {d(isub).name}';
TimeFolds(ismember(TimeFolds,{'.','..'})) = [];