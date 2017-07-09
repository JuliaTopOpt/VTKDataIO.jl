
function get_filepath_info_no_ext(filepath_no_ext)
    pattern = r"^(?P<dir>.+\\)*(?P<file>.+)$"
    _matched = match(pattern, filepath_no_ext)
    return (typeof(_matched["dir"]) == Void ? "" : _matched["dir"]), _matched["file"]
end

function write_simple_file{T<:AbstractVTKSimpleData}(vtkobject::T, filepath_or_vtmfile::Union{String, WriteVTK.MultiblockFile})
    vtkfile = nothing
    if T <: VTKUnstructuredData
        points = vtkobject.point_coords
        cells = [MeshCell(VTK_CELL_TYPE[vtkobject.cell_types[i]], vtkobject.cell_connectivity[i]) 
            for i in 1:num_of_cells(vtkobject)]
        vtkfile = vtk_grid(filepath_or_vtmfile, points, cells)
    elseif T <: VTKStructuredData
        points = vtkobject.point_coords
        vtkfile = vtk_grid(filepath_or_vtmfile, points)
    elseif T <: VTKRectilinearData
        points = vtkobject.point_coords 
        vtkfile = vtk_grid(filepath_or_vtmfile, points...)
    elseif T <: VTKUniformRectilinearData
        vtkfile = vtk_grid(filepath_or_vtmfile, vtkobject.extents..., origin=vtkobject.origin, spacing=vtkobject.spacing)
    end

    for v in keys(vtkobject.point_data)
        vtk_point_data(vtkfile, vtkobject.point_data[v], v)
    end
    for v in keys(vtkobject.cell_data)
        vtk_cell_data(vtkfile, vtkobject.cell_data[v], v)
    end

    return vtk_save(vtkfile)[1]
end

function write_blocked_file(multiblock_object::AbstractVTKMultiblockData, filename_noext::String)
    xvtm = XMLDocument()
    xroot = create_root(xvtm, "VTKFile")
    set_attribute(xroot, "type", "vtkMultiBlockDataSet")
    set_attribute(xroot, "version", "1.0")
    if ENDIAN_BOM == 0x04030201
        set_attribute(xroot, "byte_order", "LittleEndian")
    else
        set_attribute(xroot, "byte_order", "BigEndian")
    end
    new_child(xroot, "vtkMultiBlockDataSet")
    xMBDS = find_element(xroot, "vtkMultiBlockDataSet")

    filenames = []
    i = 1
    for block in simple_block_generator(multiblock_object)
        vtkFilename_noext = @sprintf("%s.z%02d", filename_noext, i)
        if isa(block, VTKPolyData)
            push!(filenames, splitdir(write_vtp(block, vtkFilename_noext))[2])            
        else
            push!(filenames, splitdir(write_simple_file(block, vtkFilename_noext))[2])
        end
        i += 1
    end

    for i in 1:length(filenames)
        fname = filenames[i]
        xBlock = new_child(xMBDS, "Block")
        set_attribute(xBlock, "index", "$i")

        xDataSet = new_child(xBlock, "DataSet")
        set_attribute(xDataSet, "index", "0")
        set_attribute(xDataSet, "file",  fname)
    end
    vtmfilepath = string(filename_noext, ".vtm")
    save_file(xvtm, vtmfilepath)

    return vtmfilepath
end

function write_timeseries_file(timeseries_object::AbstractTimeSeriesVTKData, filename_noext::String)
    xpvd = XMLDocument()
    xroot = create_root(xpvd, "VTKFile")
    set_attribute(xroot, "type", "Collection")
    set_attribute(xroot, "version", "1.0")
    if ENDIAN_BOM == 0x04030201
        set_attribute(xroot, "byte_order", "LittleEndian")
    else
        set_attribute(xroot, "byte_order", "BigEndian")
    end
    set_attribute(xroot, "compressor", "vtkZLibDataCompressor")
    new_child(xroot, "Collection")
    xMBDS = find_element(xroot, "Collection")

    main_dir, file = get_filepath_info_no_ext(filename_noext)
    subdir = main_dir*file
    mkpath(subdir)

    _length = length(timeseries_object)
    no_of_digits = 1
    while no_of_digits*10 - 1 <= _length
        no_of_digits += 1
    end

    filenames = []
    for i in 1:length(timeseries_object)
        new_filepath_no_ext = string(subdir, "\\", file, "_", "0"^(no_of_digits - length(string(i))), i)
        T = typeof(timeseries_object[i])
        if T <: VTKPolyData
            push!(filenames, write_vtp(timeseries_object[i], new_filepath_no_ext))            
        elseif T <: AbstractVTKSimpleData
            push!(filenames, write_simple_file(timeseries_object[i], new_filepath_no_ext))
        elseif T <: AbstractVTKMultiblockData
            push!(filenames, write_blocked_file(timeseries_object[i], new_filepath_no_ext))
        end
    end

    for i in 1:length(filenames)
        fname = filenames[i]
        t = timeseries_object.timemarkers[i]
        xDataSet = new_child(xMBDS, "DataSet")
        set_attribute(xDataSet, "timestep", @sprintf("%f", t))
        set_attribute(xDataSet, "part", "0")
        set_attribute(xDataSet, "file", fname)
    end

    pvdfilepath = string(filename_noext, ".pvd")
    save_file(xpvd, pvdfilepath)

    return pvdfilepath
end

function valid_to_write(vtkobject::AbstractStaticVTKData)
    T = typeof(vtkobject)
    T <: AbstractVTKMultiblockData && begin
        _out = true
        for i in vtkobject
            _T = typeof(i)
           (_out = _T.parameters[1] <: Union{Float64,Int}) || break
        end
        _out
    end || T.parameters[1] <: Union{Float64,Int} || return false, "Only float or integer data can be written."
    return true, ""
end

function valid_to_write(vtkobject::AbstractTimeSeriesVTKData)
    for b in vtkobject
       valid, _error = valid_to_write(b)
       valid ? continue : return false, _error
    end
    return true, ""
end

function write_vtk{T<:AbstractVTKData}(vtkobject::T, filepath_no_ext::String; time_resolution::Int=1, validation=false, repeat_cells=false)
    if validation
        valid, _error = is_valid(vtkobject, repeat_cells=repeat_cells)
        valid || throw("Invalid data, cannot be written: $_error")
        valid, _error = valid_to_write(vtkobject)
        valid || throw("Invalid data, cannot be written: $_error")
    end

    if T <: VTKPolyData
        return write_vtp(vtkobject, filepath_no_ext)
    elseif T <: AbstractVTKSimpleData
        return write_simple_file(vtkobject, filepath_no_ext)
    elseif T <: AbstractVTKMultiblockData
        return write_blocked_file(vtkobject, filepath_no_ext)
    elseif T <: AbstractTimeSeriesVTKData
        if time_resolution == 1
            return write_timeseries_file(vtkobject, filepath_no_ext)
        else
            _vtkobject = deepcopy(vtkobject)
            increase_resolution!(_vtkobject, time_resolution)
            return write_timeseries_file(_vtkobject, filepath_no_ext)
        end
    end
    throw("Type of object input is not supported.")
end

