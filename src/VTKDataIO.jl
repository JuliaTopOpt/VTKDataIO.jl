module VTKDataIO

using VTKDataTypes
using WriteVTK
using PyCall
using LightXML
using Iterators

@pyimport vtk.util.numpy_support as vtkns
@pyimport vtk as vtk
@pyimport numpy as np

include("vtkreaders.jl")
include("vtkwriters.jl")
include("vtpwriter.jl")
include("stlIO.jl")
include("extra_readers.jl")
include("PyVTK.jl")
include("visualize.jl")

export read_vtk, write_vtk, read_stl, write_stl, valid_to_write, PyVTK, 
    _VTKDataTypes, visualize, save_x3d, read_obj, read_ply

end
