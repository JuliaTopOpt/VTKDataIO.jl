
vtk_to_julia(x, T) = Array{T}(PyReverseDims(vtkns.vtk_to_numpy(x)))

function pytype_as_string(x) 
    pytype = pybuiltin(:str)(pybuiltin(:type)(x))
    if contains(pytype, "<type '")
        return replace(replace(pytype, "<type '", ""), "'>", "")
    elseif contains(pytype, "<class '")
        return replace(replace(pytype, "<class '", ""), "'>", "")
    end
end

function read_vtk(reader::PyCall.PyObject)
    if in(reader[:GetFileName](), readdir())
        reader[:Update]()
    else
        throw("Empty reader object or file name is not in the current directory.")
    end
    (reader[:GetClassName]() == "vtkXMLUnstructuredGridReader" ||
    reader[:GetClassName]() == "vtkXMLStructuredGridReader" ||
    reader[:GetClassName]() == "vtkXMLRectilinearGridReader" ||
    reader[:GetClassName]() == "vtkXMLImageDataReader" ||
    reader[:GetClassName]() == "vtkXMLMultiBlockDataReader" ||
    reader[:GetClassName]() == "vtkGenericDataObjectReader") && return extract_static_data(reader)

    throw("Reader type not supported yet.")
end

function read_vtk(filepath::String)
    pattern = r"\.(?P<ext>\w+)$"
    _match = match(pattern, filepath)
    ext = _match["ext"]
    if ext == "pvd"
        read_timeseries_pvd(filepath)
    #=
    elseif "ex2" == ext
        read_timeseries_ex2(filepath)
    =#
    else
        read_static_vtk(filepath)
    end
end

function read_timeseries_pvd(pvd_filepath)
    xdoc = parse_file(pvd_filepath)
    xroot = root(xdoc)
    datasets_iter = child_elements(xroot["Collection"][1])
    timemarkers = [parse(attribute(i, "timestep")) for i in datasets_iter]
    filepaths = [attribute(i, "file") for i in datasets_iter]

    static_datasets = [read_static_vtk(f) for f in filepaths]
    dataset = VTKTimeSeriesData(timemarkers, static_datasets)
end

function read_static_vtk(filepath)
    pattern = r"\.(?P<ext>\w+)$"
    _match = match(pattern, filepath)
    ext = _match["ext"]
    if "vtu" == ext
        reader = vtk.vtkXMLUnstructuredGridReader()
    elseif "vts" == ext
        reader = vtk.vtkXMLStructuredGridReader()
    elseif "vtr" == ext
        reader = vtk.vtkXMLRectilinearGridReader()        
    elseif "vti" == ext
        reader = vtk.vtkXMLImageDataReader()        
    elseif "vtm" == ext
        reader = vtk.vtkXMLMultiBlockDataReader()
    elseif "vtk" == ext
        reader = vtk.vtkGenericDataObjectReader()
    elseif "vtp" == ext
        reader = vtk.vtkXMLPolyDataReader()
#=
    elseif "e" == ext
        reader = vtk.vtkExodusReader()
    elseif "stl" == ext
        reader = vtk.vtkSTLReader()
    elseif "obj" == ext
        reader = vtk.vtkOBJReader() 
    elseif "3ds" == ext
        reader = vtk.vtk3DSImporter()
    elseif "vrml" == ext
        reader = vtk.vtkVRMLImporter()
=#
    else
        throw("Format is not supported yet.")
    end

    reader[:SetFileName](filepath)
    try
        reader[:Update]()
    catch
        throw("File corrupt.")
    end
    return extract_static_data(reader)
end

function extract_static_data(reader)
    block = nothing
    try
        block = reader[:GetOutput]()
    catch
        throw("File corrupt.")
    end
    if block[:GetClassName]() == "vtkMultiBlockDataset"
        return extract_blocked_data(block)
    else
        return extract_simple_data(block)
    end
end

function extract_simple_data(block)
    data_type = block[:GetClassName]()
    if data_type == "vtkUnstructuredGrid"
        return extract_unstructured_data(block)
    elseif data_type == "vtkPolyData"
        return extract_unstructured_data(block, true)
    elseif data_type == "vtkStructuredGrid"
        return extract_structured_data(block)
    elseif data_type == "vtkRectilinearGrid"
        return extract_rectilinear_data(block)
    elseif data_type == "vtkImageData"
        return extract_image_data(block)
    elseif data_type == "vtkMultiBlockDataSet"
        return extract_blocked_data(block)
    else
        throw("Format not supported yet.")
    end
end

function get_unstructured_point_and_cell_data(block)
    point_data = Dict{String, Array{Float64}}()
    point_vars_names = []
    for i in 1:block[:GetPointData]()[:GetNumberOfArrays]()
        var_name = block[:GetPointData]()[:GetArrayName](i-1)
        push!(point_vars_names, var_name)
        _point_data = block[:GetPointData]()[:GetArray](i-1)
        if pytype_as_string(_point_data) != "NoneType"
            point_data[var_name] = vtk_to_julia(_point_data, Float64)
        end
    end

    cell_data = Dict{String, Array{Float64}}()
    cell_vars_names = []
    for i in 1:block[:GetCellData]()[:GetNumberOfArrays]()
        var_name = block[:GetCellData]()[:GetArrayName](i-1)
        push!(cell_vars_names, var_name)
        _cell_data = block[:GetCellData]()[:GetArray](i-1)
        if pytype_as_string(_cell_data) != "NoneType"
            cell_data[var_name] = vtk_to_julia(_cell_data, Float64)
        end
    end
    return point_data, cell_data
end

function extract_unstructured_data(block, poly = false)
    point_coords = vtk_to_julia(block[:GetPoints]()[:GetData](), Float64)
    num_of_cells = block[:GetNumberOfCells]()

    _cell_types = Vector{Int}()
    if !poly
        _cell_types = vtk_to_julia(block[:GetCellTypesArray](), Int)
        vtk_cell_conn = vtk_to_julia(block[:GetCells]()[:GetData](), Int) 
        _cell_connectivity = Vector{Int}[]
        i = 1
        while i <= length(vtk_cell_conn)
            push!(_cell_connectivity, vtk_cell_conn[i+1:i+vtk_cell_conn[i]] .+ 1)
            i = i + vtk_cell_conn[i] + 1
        end
    else
        poly_pyclasses = [5,7,9]
        vert_pyclasses = [1,2]
        line_pyclasses = [3,4]
        strip_pyclasses = [6]

        _cell_types = Int[]
        _cell_connectivity = Vector{Int}[]

        vtk_verts = vtk_to_julia(block[:GetVerts]()[:GetData](), Int)
        i = 1
        while i <= length(vtk_verts)
            npoints = vtk_verts[i]
            if npoints == 1
                push!(_cell_types, 1)
            else
                push!(_cell_types, 2)
            end
            push!(_cell_connectivity, vtk_verts[i+1:i+npoints] .+ 1)
            i += npoints + 1
        end

        vtk_lines = vtk_to_julia(block[:GetLines]()[:GetData](), Int)
        i = 1
        while i <= length(vtk_lines)
            npoints = vtk_lines[i]
            if npoints == 2
                push!(_cell_types, 3)
            else
                push!(_cell_types, 4)
            end
            push!(_cell_connectivity, vtk_lines[i+1:i+npoints] .+ 1)
            i += npoints + 1
        end

        vtk_polys = vtk_to_julia(block[:GetPolys]()[:GetData](), Int)
        i = 1
        while i <= length(vtk_polys)
            npoints = vtk_polys[i]
            if npoints == 3
                push!(_cell_types, 5)
            elseif npoints == 4
                push!(_cell_types, 9)
            else
                push!(_cell_types, 7)
            end
            push!(_cell_connectivity, vtk_polys[i+1:i+npoints] .+ 1)
            i += npoints + 1
        end

        vtk_strips = vtk_to_julia(block[:GetStrips]()[:GetData](), Int)
        i = 1
        while i <= length(vtk_strips)
            npoints = vtk_strips[i]
            push!(_cell_types, 6)
            push!(_cell_connectivity, vtk_strips[i+1:i+npoints] .+ 1)
            i += npoints + 1
        end
    end

    point_data, cell_data = get_unstructured_point_and_cell_data(block)

    if poly
        return VTKPolyData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
    end

    return VTKUnstructuredData(point_coords, _cell_types, _cell_connectivity, point_data, cell_data)
end

function get_structured_point_and_cell_data(block)
    _extents = block[:GetDimensions]()
    point_data = Dict{String, Array{Float64}}()
    point_vars_names = []
    for i in 1:block[:GetPointData]()[:GetNumberOfArrays]()
        var_name = block[:GetPointData]()[:GetArrayName](i-1)
        push!(point_vars_names, var_name)
        _point_data = block[:GetPointData]()[:GetArray](i-1)
        if pytype_as_string(_point_data) != "NoneType"
            _point_data_ = vtk_to_julia(_point_data, Float64)
            var_dim = length(size(_point_data_)) == 1 ? 1 : size(_point_data_,1)
            if var_dim == 1
                point_data[var_name] = reshape(_point_data_, _extents)
            else
                point_data[var_name] = reshape(_point_data_, (var_dim, _extents...))
            end
        end
    end

    cell_extents = (([_extents...] .- 1)...)
    cell_data = Dict{String, Array{Float64}}()
    cell_vars_names = []
    for i in 1:block[:GetCellData]()[:GetNumberOfArrays]()
        var_name = block[:GetCellData]()[:GetArrayName](i-1)
        push!(cell_vars_names, var_name)
        _cell_data = block[:GetCellData]()[:GetArray](i-1)
        if pytype_as_string(_cell_data) != "NoneType"
            _cell_data_ = vtk_to_julia(_cell_data, Float64)
            var_dim = length(size(_cell_data_)) == 1 ? 1 : size(_cell_data_,1)
            if var_dim == 1
                cell_data[var_name] = reshape(_cell_data_, cell_extents)
            else
                cell_data[var_name] = reshape(_cell_data_, (var_dim, cell_extents...))
            end
        end
    end

    return point_data, cell_data
end

function extract_structured_data(block)
    _point_coords = vtk_to_julia(block[:GetPoints]()[:GetData](), Float64)

    _extents = block[:GetDimensions]()
    _dim = length(_extents)
    point_coords = reshape(_point_coords, (_dim, _extents...))
    point_data, cell_data = get_structured_point_and_cell_data(block)
    return VTKStructuredData(point_coords, point_data, cell_data)
end
function extract_rectilinear_data(block)
    extents = block[:GetDimensions]()
    dim = length(extents)
    num_of_points = prod(extents)

    point_coords = Vector{Float64}[]
    for j in 1:dim
        push!(point_coords, zeros(extents[j]))
    end
    for j in 1:dim
        k = 1
        for i in 1:prod(extents[1:j-1]):prod(extents[1:j])
            point_coords[j][k] = block[:GetPoint](i-1)[j]
            k += 1
        end
    end
    
    point_data, cell_data = get_structured_point_and_cell_data(block)
    return VTKRectilinearData(point_coords, point_data, cell_data)
end
function extract_image_data(block)
    origin = [block[:GetOrigin]()...]
    spacing = [block[:GetSpacing]()...]
    extents = [block[:GetDimensions]()...]
    point_extents = (extents...)
    cell_extents = (([extents...] .- 1)...)
    point_data, cell_data = get_structured_point_and_cell_data(block)
    return VTKUniformRectilinearData(origin, spacing, extents, point_data, cell_data)
end

function get_blocks(multiblock_pyobject)
    blocks = []
    if pytype_as_string(multiblock_pyobject) == "NoneType"
        nothing
    elseif multiblock_pyobject[:GetClassName]() != "vtkMultiBlockDataSet"
        push!(blocks, multiblock_pyobject)
    else
        for i in 1:multiblock_pyobject[:GetNumberOfBlocks]()
            for block in get_blocks(multiblock_pyobject[:GetBlock](i-1))
                if pytype_as_string(block) == "NoneType"
                    continue
                else
                    push!(blocks, block)
                end
            end
        end
    end
    return blocks
end

function extract_blocked_data(multiblock_pyobject)
    blocks = get_blocks(multiblock_pyobject)
    multiblock_data = AbstractStaticVTKData{Float64}[]
    for block in blocks
        block_data = extract_simple_data(block)
        push!(multiblock_data, block_data)
    end
    return VTKMultiblockData(multiblock_data)
end

