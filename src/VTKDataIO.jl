module VTKDataIO

__precompile__(false)

using VTKDataTypes
using WriteVTK
using PyCall
using LightXML
using IterTools
using Printf

if get(ENV, "JULIA_REGISTRYCI_AUTOMERGE", "false") == "false"
    @eval begin
        const vtkns = pyimport_conda("vtk.util.numpy_support", "vtk")
        const vtk = pyimport_conda("vtk", "vtk")
        const np = pyimport_conda("numpy", "numpy")

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
end

end
