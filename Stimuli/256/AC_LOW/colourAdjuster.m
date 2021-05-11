list = dir('*.jpeg');

for i = 1:length(list)
    oldName = list(i,1).name;
    img = imread(oldName);
    imgCol = colourblind(img, 'protanopic');
    newName = replace(oldName, '.jpeg','_col.jpeg')    
    imwrite(imgCol,newName);
end
