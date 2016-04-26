function outVec = parseXML(fileName,selectedOutput)

try
   tree = xmlread(fileName);
catch
   error('Failed to read XML file %s.',filename);
end

if isempty(selectedOutput)
    outVec = [];
    return;
end

numSteps = tree.getElementsByTagName('d-0-0').getLength;
outVec = zeros(size(selectedOutput,1),numSteps);

for i=0:numSteps-1
    valueNum = strcat('value',num2str(i));
    
    valueChildren = tree.getElementsByTagName(valueNum).item(0).getChildNodes;
    
    for j=1:size(selectedOutput,1)
        dim = ['d-' num2str(selectedOutput(j,1)) '-' num2str(selectedOutput(j,2))];
        outVec(j,i+1) = str2double(valueChildren.getElementsByTagName(dim).item(0).getChildNodes.item(0).getData);
    end
    
end

end
