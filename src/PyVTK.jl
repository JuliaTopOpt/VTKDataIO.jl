VTKClasses = Dict(0=>vtk.vtkEmptyCell,
        1=>vtk.vtkVertex, 
        2=>vtk.vtkPolyVertex,
        3=>vtk.vtkLine, 
        4=>vtk.vtkPolyLine,
        5=>vtk.vtkTriangle,
        6=>vtk.vtkTriangleStrip,
        7=>vtk.vtkPolygon,
        8=>vtk.vtkPixel,
        9=>vtk.vtkQuad,
        10=>vtk.vtkTetra,
        11=>vtk.vtkVoxel,
        12=>vtk.vtkHexahedron,
        13=>vtk.vtkWedge,
        14=>vtk.vtkPyramid,
        15=>vtk.vtkPentagonalPrism,
        16=>vtk.vtkHexagonalPrism,
        21=>vtk.vtkQuadraticEdge, 
        22=>vtk.vtkQuadraticTriangle, 
        23=>vtk.vtkQuadraticQuad,
        24=>vtk.vtkQuadraticTetra,
        25=>vtk.vtkQuadraticHexahedron,
        26=>vtk.vtkQuadraticWedge,
        27=>vtk.vtkQuadraticPyramid,
        28=>vtk.vtkBiQuadraticQuad,
        29=>vtk.vtkTriQuadraticHexahedron,
        30=>vtk.vtkQuadraticLinearQuad,
        31=>vtk.vtkQuadraticLinearWedge,
        32=>vtk.vtkBiQuadraticQuadraticWedge,
        33=>vtk.vtkBiQuadraticQuadraticHexahedron,
        34=>vtk.vtkBiQuadraticTriangle,
        35=>vtk.vtkCubicLine,
        41=>vtk.vtkConvexPointSet,
        42=>vtk.vtkPolyhedron)

function PyVTK(dataset::VTKUnstructuredData)
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    if dim(dataset) == 2
        points[:SetData](vtkns.numpy_to_vtk(PyReverseDims([dataset.point_coords; zeros(num_of_points(dataset))'])))
    else
        points[:SetData](vtkns.numpy_to_vtk(PyReverseDims(dataset.point_coords)))
    end
    grid[:SetPoints](points)

    jl_cells = Int[]
    jl_cell_locs = Int[]
    k = 0
    for i in 1:num_of_cells(dataset)
        push!(jl_cell_locs, k)
        push!(jl_cells, length(dataset.cell_connectivity[i]))
        k += 1

        if dataset.cell_types[i] == 8
            pixel_cc = dataset.cell_connectivity[i]
            push!(jl_cells, pixel_cc[1]-1)
            push!(jl_cells, pixel_cc[2]-1)
            push!(jl_cells, pixel_cc[4]-1)
            push!(jl_cells, pixel_cc[3]-1)
            k += 4
        elseif dataset.cell_types[i] == 11
            voxel_cc = dataset.cell_connectivity[i]
            push!(jl_cells, voxel_cc[1]-1)
            push!(jl_cells, voxel_cc[2]-1)
            push!(jl_cells, voxel_cc[4]-1)
            push!(jl_cells, voxel_cc[3]-1)
            push!(jl_cells, voxel_cc[5]-1)
            push!(jl_cells, voxel_cc[6]-1)
            push!(jl_cells, voxel_cc[7]-1)
            push!(jl_cells, voxel_cc[8]-1)
            k += 8
        else
            for j in dataset.cell_connectivity[i]
                push!(jl_cells, j-1)
                k += 1
            end            
        end
    end
    _cells = vtkns.numpy_to_vtk(np.array(jl_cells), 0, 12) #vtkTypeIdArray() of serialized cell connectivity
    ncells = num_of_cells(dataset)
    vtkCells = vtk.vtkCellArray()
    vtkCells[:SetCells](ncells, _cells)
    vtk_cell_types = vtkns.numpy_to_vtk(np.array(UInt8.(dataset.cell_types)), 0, 3) #vtkUnsignedCharArray() of cell types 
    vtk_cell_locs = vtkns.numpy_to_vtk(np.array(jl_cell_locs), 0, 12) #vtkTypeIdArray() of cell locations
    grid[:SetCells](vtk_cell_types, vtk_cell_locs, vtkCells)

    for m in keys(dataset.point_data)
        vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(dataset.point_data[m]))
        vtk_data_array[:SetName](m)
        grid[:GetPointData]()[:AddArray](vtk_data_array)
    end

    for m in keys(dataset.cell_data)
        vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(dataset.cell_data[m]))
        vtk_data_array[:SetName](m)
        grid[:GetCellData]()[:AddArray](vtk_data_array)
    end

    return grid
end

function PyVTK(dataset::VTKPolyData)
    poly_pyclasses = [5,7,8,9]
    vert_pyclasses = [1,2]
    line_pyclasses = [3,4]
    strip_pyclasses = [6]

    new_inds = Int[]
    polydata = vtk.vtkPolyData()

    points = vtk.vtkPoints()
    if dim(dataset) == 2
        points[:SetData](vtkns.numpy_to_vtk(PyReverseDims([dataset.point_coords; zeros(num_of_points(dataset))'])))
    else
        points[:SetData](vtkns.numpy_to_vtk(PyReverseDims(dataset.point_coords)))
    end        
    polydata[:SetPoints](points)

    jl_verts = Int[]
    jl_lines = Int[]
    jl_polys = Int[]
    jl_strips = Int[]
    vertinds = Int[]
    lineinds = Int[]
    polyinds = Int[]
    stripinds = Int[]

    for i in 1:num_of_cells(dataset)
        _cell_type = dataset.cell_types[i]
        if _cell_type ∈ vert_pyclasses
            push!(vertinds, i)
            push!(jl_verts, length(dataset.cell_connectivity[i]))
            for j in dataset.cell_connectivity[i]
                push!(jl_verts, j-1)
            end
        elseif _cell_type ∈ line_pyclasses
            push!(lineinds, i)
            push!(jl_lines, length(dataset.cell_connectivity[i]))
            for j in dataset.cell_connectivity[i]
                push!(jl_lines, j-1)
            end            
        elseif _cell_type ∈ poly_pyclasses
            push!(polyinds, i)
            push!(jl_polys, length(dataset.cell_connectivity[i]))
            if _cell_type == 8
                pixel_cc = dataset.cell_connectivity[i]
                push!(jl_polys, pixel_cc[1]-1)
                push!(jl_polys, pixel_cc[2]-1)
                push!(jl_polys, pixel_cc[4]-1)
                push!(jl_polys, pixel_cc[3]-1) 
            else
                for j in dataset.cell_connectivity[i]
                    push!(jl_polys, j-1)
                end        
            end
        elseif _cell_type ∈ strip_pyclasses
            push!(stripinds, i)
            push!(jl_strips, length(dataset.cell_connectivity[i]))
            for j in dataset.cell_connectivity[i]
                push!(jl_strips, j-1)
            end
        else
            throw("Cell $i is not supported.")
        end
    end

    _verts = vtkns.numpy_to_vtk(np.array(jl_verts), 0, 12) #vtkTypeIdArray() of serialized cell connectivity
    nverts = length(vertinds)
    verts = vtk.vtkCellArray()
    verts[:SetCells](nverts, _verts)
    polydata[:SetVerts](verts)

    _lines = vtkns.numpy_to_vtk(np.array(jl_lines), 0, 12) #vtkTypeIdArray() of serialized cell connectivity
    nlines = length(lineinds)
    lines = vtk.vtkCellArray()
    lines[:SetCells](nlines, _lines)
    polydata[:SetLines](lines)

    _polys = vtkns.numpy_to_vtk(np.array(jl_polys), 0, 12) #vtkTypeIdArray() of serialized cell connectivity
    npolys = length(polyinds)
    polys = vtk.vtkCellArray()
    polys[:SetCells](npolys, _polys)
    polydata[:SetPolys](polys)

    _strips = vtkns.numpy_to_vtk(np.array(jl_strips), 0, 12) #vtkTypeIdArray() of serialized cell connectivity
    nstrips = length(stripinds)
    strips = vtk.vtkCellArray()
    strips[:SetCells](nstrips, _strips)
    polydata[:SetStrips](strips)

    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(dataset.point_data[m]))
        else
            vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(dataset.point_data[m]))
        end
        vtk_data_array[:SetName](m)
        polydata[:GetPointData]()[:AddArray](vtk_data_array)
    end

    new_cell_inds = chain(vertinds, lineinds, polyinds, stripinds)
    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            new_cell_data = zeros(num_of_cells(dataset))
            i = 1
            for j in new_cell_inds
                new_cell_data[i] = dataset.cell_data[m][j]
                i += 1
            end
            vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(new_cell_data))
        else
            _size = size(dataset.cell_data[m])
            new_cell_data = zeros(_size...)
            i = 1
            for j in new_cell_inds
                new_cell_data[:,i] = dataset.cell_data[m][:,j]
                i += 1
            end
            vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(new_cell_data))
        end            
        vtk_data_array[:SetName](m)
        polydata[:GetCellData]()[:AddArray](vtk_data_array)
    end

    return polydata
end

function PyVTK(dataset::VTKStructuredData)
    grid = vtk.vtkStructuredGrid()
    points = vtk.vtkPoints()

    _point_coords = reshape(dataset.point_coords, (dim(dataset), num_of_points(dataset)))
    if dim(dataset) == 2
        points[:SetData](vtkns.numpy_to_vtk(PyReverseDims([_point_coords; zeros(num_of_points(dataset))'])))
    else
        points[:SetData](vtkns.numpy_to_vtk(PyReverseDims(_point_coords)))
    end
    grid[:SetPoints](points)

    if _dim == 2
        grid[:SetDimensions](_point_extents..., 1)
    else
        grid[:SetDimensions](_point_extents...)
    end
    
    add_structured_data_to_grid!(grid, dataset)
    return grid
end

function PyVTK(dataset::VTKRectilinearData)
    grid = vtk.vtkRectilinearGrid()

    _dim = dim(dataset)
    _point_extents = extents(dataset)

    xCoords = vtkns.numpy_to_vtk(PyReverseDims(dataset.point_coords[1]))
    yCoords = vtkns.numpy_to_vtk(PyReverseDims(dataset.point_coords[2]))
    if _dim == 2
        zCoords = vtkns.numpy_to_vtk(PyReverseDims([0.]))
    else
        zCoords = vtkns.numpy_to_vtk(PyReverseDims(dataset.point_coords[3]))
    end

    if _dim == 2
        grid[:SetDimensions](_point_extents..., 1)
    else
        grid[:SetDimensions](_point_extents...)
    end

    grid[:SetXCoordinates](xCoords)
    grid[:SetYCoordinates](yCoords)
    grid[:SetZCoordinates](zCoords)

    add_structured_data_to_grid!(grid, dataset)
    return grid
end

function PyVTK(dataset::VTKUniformRectilinearData)
    grid = vtk.vtkImageData()

    _dim = dim(dataset)
    _point_extents = extents(dataset)

    if _dim == 2
        grid[:SetOrigin](dataset.origin..., 0)
        grid[:SetSpacing](dataset.spacing..., 0)
        grid[:SetDimensions](_point_extents..., 1)
    else
        grid[:SetOrigin](dataset.origin...)
        grid[:SetSpacing](dataset.spacing...)
        grid[:SetDimensions](_point_extents...)        
    end

    add_structured_data_to_grid!(grid, dataset)
    return grid
end

function add_structured_data_to_grid!(grid::PyObject, dataset::AbstractVTKStructuredData)
    for m in keys(dataset.point_data)
        _var_dim = var_dim(dataset, m, "Point")
        if _var_dim == 1
            data_array = reshape(dataset.point_data[m], (1, num_of_points(dataset)))
        else
            data_array = reshape(dataset.point_data[m], (_var_dim, num_of_points(dataset)))
        end
        vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(data_array))
        vtk_data_array[:SetName](m)
        grid[:GetPointData]()[:AddArray](vtk_data_array)
    end

    for m in keys(dataset.cell_data)
        _var_dim = var_dim(dataset, m, "Cell")
        if _var_dim == 1
            data_array = reshape(dataset.cell_data[m], (1, num_of_cells(dataset)))
        else
            data_array = reshape(dataset.cell_data[m], (_var_dim, num_of_cells(dataset)))
        end
        vtk_data_array = vtkns.numpy_to_vtk(PyReverseDims(data_array))
        vtk_data_array[:SetName](m)
        grid[:GetCellData]()[:AddArray](vtk_data_array)
    end
end

function PyVTK(dataset::VTKMultiblockData)
    mb = vtk.vtkMultiBlockDataSet()
    for i in 1:length(dataset)
        block = dataset[i]
        mb[:SetBlock](i-1, PyVTK(block))
    end
    return mb
end

function _VTKDataTypes(pyvtk::PyObject)
    if pyvtk[:GetClassName]() == "vtkMultiBlockDataSet"
        return extract_blocked_data(pyvtk)
    else
        return extract_simple_data(pyvtk)
    end
end

function write_vtk(pyvtk::PyObject, filepath_noext::String)
    classname = pyvtk[:GetClassName]()
    if classname == "vtkUnstructuredGrid"
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer[:SetFileName](string(filepath_noext, ".vtu"))
    elseif classname == "vtkStructuredGrid"
        writer = vtk.vtkXMLStructuredGridWriter()
        writer[:SetFileName](string(filepath_noext, ".vts"))
    elseif classname == "vtkRectilinearGrid"
        writer = vtk.vtkXMLRectilinearGridWriter()
        writer[:SetFileName](string(filepath_noext, ".vtr"))
    elseif classname == "vtkImageData"
        writer = vtk.vtkXMLImageDataWriter()
        writer[:SetFileName](string(filepath_noext, ".vti"))
    elseif classname == "vtkMultiBlockDataSet"
        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer[:SetFileName](string(filepath_noext, ".vtm"))
    else
        throw("Type of VTK object is not supported.")
    end

    writer[:SetInputData](pyvtk)
    writer[:Update]()

    return filepath_noext
end
