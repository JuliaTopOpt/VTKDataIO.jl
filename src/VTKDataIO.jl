module VTKDataIO

# __precompile__(false)

using VTKDataTypes
using WriteVTK
using PythonCall
using LightXML
using IterTools
using Printf

const vtk = PythonCall.pynew() # initially NULL
const vtkns = PythonCall.pynew() # initially NULL
const np = PythonCall.pynew() # initially NULL
function __init__()
  PythonCall.pycopy!(vtk, pyimport("vtk"))
  PythonCall.pycopy!(vtkns, pyimport("vtk.util.numpy_support"))
  PythonCall.pycopy!(np, pyimport("numpy"))
end

include("vtkreaders.jl")
include("vtkwriters.jl")
include("vtpwriter.jl")
include("stlIO.jl")
include("extra_readers.jl")
include("PyVTK.jl")
include("visualize.jl")

export read_vtk, write_vtk, read_stl, write_stl, valid_to_write, PyVTK, 
_VTKDataTypes, visualize, write_x3d, write_ply, read_obj, read_ply, visualize_3ds

end
