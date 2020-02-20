function [vertex, faces] = readAneurysmMeshObj(fileName)

OBJ=read_wobj(fileName);
faces = [];
vertex = [];

for i = 1 : numel(OBJ.objects)
    if isfield(OBJ.objects(1,i).data,'vertices')
        if numel(OBJ.objects(1,i).data.vertices) > numel(faces) 
            faces = OBJ.objects(1,i).data.vertices;
        end
    end
end

if isfield(OBJ,'vertices')
    vertex = OBJ.vertices;
end

if size(vertex,1) == 3
    vertex = vertex';
end
if size(faces,1) == 3
    faces = faces';
end