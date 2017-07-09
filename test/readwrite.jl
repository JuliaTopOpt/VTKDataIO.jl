using VTKDataTypes
using VTKDataIO

function write_read_3d_image()
    origin = [0, 0, 0]
    spacing = [1., 2., 1.,]
    _extents = [2, 3, 5]

    image = VTKImageData(origin, spacing, _extents)
    image.point_data["Point scalar data"] = rand(extents(image)...)
    image.cell_data["Cell scalar data"] = rand(cell_extents(image)...)
    image.point_data["Point vector data"] = rand(3, extents(image)...)
    image.cell_data["Cell vector data"] = rand(3, cell_extents(image)...)

    write_vtk(image, "test_3d_image") # Without extension

    _image = read_vtk("test_3d_image.vti")
    @test image == _image

    return image
end

function write_read_3d_rectilinear()
    image = write_read_3d_image()
    rectilinear = VTKRectilinearData(image)
    write_vtk(rectilinear, "test_3d_rect") # Without extension
    _rectilinear = read_vtk("test_3d_rect.vtr")
    @test rectilinear == _rectilinear
    return rectilinear
end

function write_read_3d_structured()
    rectilinear = write_read_3d_rectilinear()
    structured = VTKStructuredData(rectilinear)
    write_vtk(structured, "test_3d_struct") # Without extension
    _structured = read_vtk("test_3d_struct.vts")
    @test structured == _structured
    return structured
end

function write_read_3d_unstructured()
    structured = write_read_3d_structured()
    unstruct = VTKUnstructuredData(structured)
    write_vtk(unstruct, "test_3d_unstruct") # Without extension
    _unstruct = read_vtk("test_3d_unstruct.vtu")
    @test unstruct == _unstruct
    return unstruct
end

function write_read_3d_polydata()
    rectilinear = write_read_3d_rectilinear()
    polydata = VTKPolyData(rectilinear)
    write_vtk(polydata, "test_3d_poly") # Without extension
    #When writing 2D polydata, it is automatically changed to 3D
    _polydata = read_vtk("test_3d_poly.vtp")
    @test polydata == _polydata
    return polydata
end

function write_read_3d_multiblock()
    mb = VTKMultiblockData([write_read_3d_image(), write_read_3d_rectilinear(), 
        write_read_3d_structured(), write_read_3d_unstructured(),
        write_read_3d_polydata()]);
    write_vtk(mb, "test_3d_multiblock") # Without extension
    _mb = read_vtk("test_3d_multiblock.vtm")
    @test mb == _mb
    return mb
end

function write_read_3d_timeseries()
    #Making a timeseries data

    x = y = z = [0., 1., 2.];
    point_coords = [x,y,z];
    point_data = Dict{String, Array{Float64}}();
    cell_data = Dict{String, Array{Float64}}();
    a = VTKRectilinearData(point_coords, point_data, cell_data);
    a.cell_data["Color"] = reshape(rand(num_of_cells(a)), cell_extents(a));
    b = VTKPolyData(a);

    x = y = z = [3., 4., 5.];
    point_coords = [x,y,z];
    c = VTKPolyData(VTKRectilinearData(point_coords, point_data, cell_data));
    c.cell_data["Color"] = rand(num_of_cells(c));

    #Multiblock data
    m = VTKMultiblockData([b, c]);
    timeseries = VTKTimeSeriesData([0.], [m]);

    # Adding new time steps
    no_of_timesteps = 5;
    timestep = 0.5;
    speed = 1;
    i = 0
    while i < no_of_timesteps
        i += 1;
        new_data = deepcopy(timeseries.data[i]);

        #Random walk of first dataset
        new_data[1].point_coords = new_data[1].point_coords + 
            vcat(fill(speed*rand()-speed/2, (1,num_of_points(new_data[1]))),  
            fill(speed*rand()-speed/2, (1,num_of_points(new_data[1]))), 
            fill(speed*rand()-speed/2, (1,num_of_points(new_data[1]))));
                
        #Random walk of second dataset
        new_data[2].point_coords = new_data[2].point_coords + 
            vcat(fill(speed*rand()-speed/2, (1,num_of_points(new_data[2]))), 
            fill(speed*rand()-speed/2, (1,num_of_points(new_data[2]))), 
            fill(speed*rand()-speed/2, (1,num_of_points(new_data[2]))));

        insert_timed_data!(timeseries, i*timestep, new_data);
    end

    #Writing and reading

    write_vtk(timeseries, "test_3d_timeseries") # Without extension
    _timeseries = read_vtk("test_3d_timeseries.pvd")
    @test timeseries == _timeseries
    return timeseries
end

write_read_3d_image();
write_read_3d_rectilinear();
write_read_3d_structured();
write_read_3d_unstructured();
write_read_3d_polydata();
write_read_3d_multiblock();
write_read_3d_timeseries();
