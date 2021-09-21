function read_obj(filepath::String)
    if splitext(filepath)[2] == ".obj"
        reader = vtk.vtkOBJReader()
        reader.SetFileName(filepath)
        reader.Update()

        polydata = reader.GetOutput()
        return extract_simple_data(polydata)
    else
        throw("The file input is not an obj file.")
    end
end

function read_ply(filepath::String)
    if splitext(filepath)[2] == ".ply"
        reader = vtk.vtkPLYReader()
        reader.SetFileName(filepath)
        reader.Update()

        polydata = reader.GetOutput()
        return extract_simple_data(polydata)
    else
        throw("The file input is not a ply file.")
    end
end
