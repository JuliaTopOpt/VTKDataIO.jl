
function write_vtp(dataset::VTKPolyData, filepath_no_ext::String)
    polydata = PyVTK(dataset)
    return write_vtp(polydata, filepath_no_ext)
end

function write_vtp(polydata::PyObject, filepath_no_ext::String)
    writer = vtk.vtkXMLPolyDataWriter()
    writer[:SetFileName](string(filepath_no_ext, ".vtp"))
    writer[:SetDataModeToBinary]()
    writer[:SetDataModeToAppended]()
    writer[:SetCompressorTypeToZLib]()
    writer[:SetInputData](polydata)
    writer[:Update]()
    return string(filepath_no_ext, ".vtp")
end

