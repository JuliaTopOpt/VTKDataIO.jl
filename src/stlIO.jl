
function write_stl(dataset::AbstractStaticVTKData, filepath_no_ext::String)
    polydata = PyVTK(trianulate(dataset))
    return write_stl(polydata, filepath_no_ext)
end

function write_stl(polydata::PyCall.PyObject, filepath_no_ext::String)
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(string(filepath_no_ext, ".stl"))
    writer.SetInputData(polydata)
    writer.Update()
    return string(filepath_no_ext, ".stl")
end

function read_stl(filepath::String)
    if splitext(filepath)[2] == ".stl"
        reader = vtk.vtkSTLReader()
        reader.SetFileName(filepath)
        reader.Update()

        polydata = reader.GetOutput()
        return extract_simple_block(polydata)
    else
        throw("The file input is not an stl file.")
    end
end
