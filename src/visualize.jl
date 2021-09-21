
function visualize(dataset::AbstractStaticVTKData; color="", RGB=false, component=-1, opacity=1., window_size=(1000,600), ncolors=100, legend=false, legendtitle=color, scale_factor=1, representation=:simple)
    renderwindow, interactor = make_render_window_and_interactor(dataset, color, RGB, component, opacity, window_size, ncolors, legend, legendtitle, scale_factor, representation)
    renderwindow.Render()
    interactor.Start()
    renderwindow.Finalize()
end

function visualize_3ds(filename::AbstractString; window_size=(1000,600), background=(0.1, 0.2, 0.4), camera_position=(0, 1, 0), camera_focal=(0, 0, 0), camera_up=(0,0,1), camera_dolly=1.4)
    importer = vtk.vtk3DSImporter()
    importer.ComputeNormalsOn()
    importer.SetFileName(filename)
    importer.Read()
    renderwindow = importer.GetRenderWindow()
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderwindow)
    renderwindow.SetSize(window_size...);

    renderer = importer.GetRenderer()
    renderer.SetBackground(background);

    camera = renderer.GetActiveCamera()
    camera.SetFocalPoint(camera_focal) #Center
    camera.Azimuth(37.5)
    camera.Elevation(30)
    camera.Dolly(camera_dolly)

    renderer.ResetCamera();
    renderer.ResetCameraClippingRange();

    axes = vtk.vtkAxesActor() 
    widget = vtk.vtkOrientationMarkerWidget()
    widget.SetOutlineColor(0.9300,0.5700,0.1300)
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(interactor)
    widget.SetViewport(0., 0., 0.2, 0.2)
    widget.SetEnabled(1)
    widget.InteractiveOn()

    renderwindow.Render()
    interactor.Start();
end

function write_x3d(dataset::AbstractStaticVTKData, filepath_noext::String; color="", RGB=false, component=-1, opacity=1., window_size=(1000,600), ncolors=100, background=(.1, .2, .3), legend=false, legendtitle=color, scale_factor=1, representation=:simple)
    polydata = VTKPolyData(dataset)
    if color != "" && legend
        save_legend(polydata, string(filepath_noext, "_legend.png"), color=color, component=component, window_size=window_size, ncolors=ncolors, opacity=opacity, background=background, legendtitle=legendtitle)
    end

    renderwindow, interactor = make_render_window_and_interactor(polydata, color, RGB, component, opacity, window_size, ncolors, legend, legendtitle, scale_factor, representation)
    renderwindow.OffScreenRenderingOn()
    renderwindow.Render()

    writer = vtk.vtkX3DExporter()
    writer.SetFileName(string(filepath_noext, ".x3d"))
    writer.SetRenderWindow(renderwindow)
    writer.Update()

    renderwindow.Finalize()
end

function write_ply(dataset, filepath_noext; color="", component=-1, opacity=1., window_size=(800,600), ncolors=100, background=(.1, .2, .3), legend=false, legendtitle=color, scale_factor=1, representation=:simple)
    if representation == :wireframe
        throw("write_ply does not support wireframe representation.")
    end

    if color in keys(dataset.point_data)
        point_color = true
    elseif color in keys(dataset.cell_data)
        point_color = false
    end

    polydata = VTKPolyData(dataset)
    vtkobject = PyVTK(polydata)

    if color != ""
        if legend
            save_legend(polydata, string(filepath_noext, "_legend.png"), color=color, component=component, window_size=window_size, ncolors=ncolors, opacity=opacity, background=background, legendtitle=legendtitle)
        end
        jl_mapped_colors = get_jl_mapped_colors(vtkobject, polydata, color, component, opacity)[1]
    else
        jl_mapped_colors = fill(UInt8(255), (4, num_of_cells(polydata)))
    end

    if point_color
        mapped_colors = julia_to_vtk(jl_mapped_colors)
        mapped_colors.SetName("**_$(color)_**")
        mapped_colors.SetNumberOfComponents(size(jl_mapped_colors, 1))
        vtkobject.GetPointData().SetScalars(mapped_colors)
    else
        mapped_colors = julia_to_vtk(jl_mapped_colors)
        mapped_colors.SetName("**_$(color)_**")
        mapped_colors.SetNumberOfComponents(size(jl_mapped_colors, 1))
        vtkobject.GetCellData().SetScalars(mapped_colors)
    end

    if representation == :glyph
        vtkobject.GetPointData().SetActiveVectors(color)
        arrow_source = vtk.vtkArrowSource()
        glyph_filter = vtk.vtkGlyph3D()
        glyph_filter.SetSourceConnection(arrow_source.GetOutputPort())
        glyph_filter.SetScaleFactor(0.005*scale_factor)
        glyph_filter.SetColorModeToColorByScalar()
        glyph_filter.SetInputData(vtkobject)
        glyph_filter.OrientOn()
        glyph_filter.SetVectorModeToUseVector()
        glyph_filter.Update()
        writer = vtk.vtkPLYWriter()
        writer.SetFileName(string(filepath_noext,".ply"))
        writer.SetInputData(glyph_filter.GetOutput())
        writer.SetArrayName("**_$(color)_**")
        writer.Update()
    else
        writer = vtk.vtkPLYWriter()
        writer.SetFileName(string(filepath_noext,".ply"))
        writer.SetInputData(vtkobject)
        writer.SetArrayName("**_$(color)_**")
        writer.Update()
    end
end

function get_jl_mapped_colors(vtkobject, dataset::T, color_variable_name, component, opacity) where {T<:AbstractVTKSimpleData}
    if color_variable_name in keys(dataset.point_data)
        point_color = true
        _size = num_of_points(dataset)
        _color_variable = dataset.point_data[color_variable_name]
        _var_dim = var_dim(dataset, color_variable_name, "Point")
    elseif color_variable_name in keys(dataset.cell_data)
        point_color = false
        _size = num_of_cells(dataset)
        _color_variable = dataset.cell_data[color_variable_name]
        _var_dim = var_dim(dataset, color_variable_name, "Cell")
    end

    if T <: AbstractVTKStructuredData
        if _var_dim == 1
            color_variable = reshape(_color_variable, (_size,))
        else
            color_variable = reshape(_color_variable, (_var_dim, _size))
        end
    else
        color_variable = _color_variable
    end

    if component == -1 || component == 1 && _var_dim == 1
        if _var_dim == 1
            _values = color_variable
            cmin = min(_values...)
            cmax = max(_values...)
        else
            _values = [norm(color_variable[:,i]) for i in 1:_size]
            cmin = min(_values...)
            cmax = max(_values...)
        end
    else
        if 1 <= component <= _var_dim
            _values = color_variable[component,:]
            cmin = min(_values...)
            cmax = max(_values...)
        else
            throw("Cannot use component $component of a $color. $color only has $_var_dim components.")
        end
    end
    jl_mapped_colors = zeros(UInt8, 4, _size)
    if !(cmax ≈ cmin)
        scaled_values = (_values .- cmin) ./ (cmax - cmin)
        function fill_jl_mapped_colors(i)
            jl_mapped_colors[:,i] = UInt8[round.(color_map(scaled_values[i]) .* 255) ; round(255*opacity)]
        end
        map(fill_jl_mapped_colors, 1:_size)
    end
    return jl_mapped_colors, cmin, cmax
end

function save_legend(dataset::AbstractStaticVTKData, filepath_noext::String; color="", component=-1, window_size=(1000,600), ncolors=100, opacity=1., background=(.1, .2, .3), legendtitle=color)
    if color != ""
        renderwindow = make_legend_render_window(dataset, color, background, ncolors, legendtitle, component, opacity, window_size)
        imagefilter = vtk.vtkWindowToImageFilter()
        imagefilter.SetInput(renderwindow)
        imagefilter.SetViewport(0.8, 0, 1., 1.)

        imagewriter = vtk.vtkPNGWriter()
        imagewriter.SetFileName(string(filepath_noext,".png"))
        imagewriter.SetInputConnection(imagefilter.GetOutputPort())
        imagewriter.Write()

        renderwindow.Finalize()
    end
end

function make_legend_render_window(dataset, color, background, ncolors, legendtitle, component, opacity, window_size)
    if color in keys(dataset.point_data)
        point_color = true
        _size = num_of_points(dataset)
        _color_variable = dataset.point_data[color] 
        _var_dim = var_dim(dataset, color, "Point")
    elseif color in keys(dataset.cell_data)
        point_color = false
        _size = num_of_cells(dataset)
        _color_variable = dataset.cell_data[color] 
        _var_dim = var_dim(dataset, color, "Cell")
    end
    if component == -1 || component == 1 && _var_dim == 1
        if _var_dim == 1
            _values = reshape(_color_variable, (_size))
            cmin = min(_values...)
            cmax = max(_values...)
        else
            color_variable = reshape(_color_variable, (_var_dim, _size))
            _values = [norm(color_variable[:,i]) for i in 1:_size]
            cmin = min(_values...)
            cmax = max(_values...)
        end
    else
        if 1 <= component <= _var_dim
            color_variable = reshape(_color_variable, (_var_dim, _size))
            _values = color_variable[component,:]
            cmin = min(_values...)
            cmax = max(_values...)
        else
            throw("Cannot use component $component of a $color. $color only has $_var_dim components.")
        end
    end

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(background...)
    renderwindow = vtk.vtkRenderWindow()
    renderwindow.SetSize(window_size...)
    renderwindow.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderwindow)

    add_scalar_bar!(interactor, cmin, cmax, ncolors, opacity, legendtitle)

    renderwindow.OffScreenRenderingOn()
    renderwindow.Render()

    return renderwindow
end

function add_scalar_bar!(interactor, cmin, cmax, ncolors, opacity, legendtitle)
    jl_legend_colors = zeros(UInt8, 4, ncolors)
    function fill_jl_legend_colors(i)
        jl_legend_colors[:,i] = UInt8[round.(color_map(0+(i-1)*1/ncolors).*255); round(255*opacity)]
    end
    map(fill_jl_legend_colors, 1:ncolors)
    legend_colors = julia_to_vtk(jl_legend_colors)
 
    lut = vtk.vtkLookupTable()
    lut.SetTableRange(cmin, cmax)
    lut.SetNumberOfColors(ncolors)
    lut.SetTable(legend_colors)

    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetNumberOfLabels(8)
    scalar_bar.SetLabelFormat("%+#6.2e")
    scalar_bar.SetLookupTable(lut)
    scalar_bar.GetLabelTextProperty().SetFontFamilyToCourier()
    scalar_bar.GetLabelTextProperty().SetJustificationToRight()
    scalar_bar.GetLabelTextProperty().SetVerticalJustificationToBottom()
    scalar_bar.GetLabelTextProperty().BoldOff()
    scalar_bar.GetLabelTextProperty().ItalicOff()
    scalar_bar.GetLabelTextProperty().ShadowOff()
    scalar_bar.GetLabelTextProperty().SetColor(1, 1, 1)
    scalar_bar.SetTitle(legendtitle)
    scalar_bar.GetTitleTextProperty().SetFontFamilyToArial()
    scalar_bar.GetTitleTextProperty().ItalicOff()
    scalar_bar.GetTitleTextProperty().ShadowOff()

    # create the scalar bar widget
    scalar_bar_widget = vtk.vtkScalarBarWidget()
    scalar_bar_widget.SetInteractor(interactor)
    scalar_bar_widget.SetScalarBarActor(scalar_bar)
    scalar_bar_widget.On()
end

function color_map(v)
    if 0.0 <= v < 0.25
        b = 1.0
        g = v/0.25
        r = 0.0
    elseif 0.25 <= v < 0.5
        b = 1.0-(v-0.25)/0.25
        g = 1.0
        r = 0.0
    elseif 0.5 <= v < 0.75
        b = 0.0
        g = 1.0
        r = (v-0.5)/0.25
    else
        b = 0.0
        g = 1.0-(v-0.75)/0.25
        r = 1.0
    end
    return [r,g,b]
end

function setup_camera(camera, dataset)
    _bb = bb(dataset)
    calib = max([(_bb[i+1]-_bb[i]) for i in 1:2:length(_bb)]...)
    camera.SetFocalPoint(pseudo_center(dataset)...) #Center
    camera.Azimuth(37.5)
    camera.Elevation(30)
    camera.Dolly(0.2/calib)
end

function make_render_window_and_interactor(dataset::T, color, 
    RGB, component, opacity, window_size, ncolors, legend, legendtitle, scale_factor, 
    representation) where {T<:AbstractStaticVTKData}

    vtkobject = PyVTK(dataset)
    if color in keys(dataset.point_data)
        if RGB
            jl_mapped_colors = UInt8.(round.(reshape(dataset.point_data[color], 
                (size(dataset.point_data[color], 1), num_of_points(dataset)))))
            mapped_colors = julia_to_vtk(jl_mapped_colors)
            mapped_colors.SetNumberOfComponents(size(dataset.point_data[color], 1))
            vtkobject.GetPointData().SetScalars(mapped_colors)
        else
            jl_mapped_colors, cmin, cmax = get_jl_mapped_colors(vtkobject, dataset, color, component, opacity)
            if !(cmax ≈ cmin)
                mapped_colors = julia_to_vtk(jl_mapped_colors)
                vtkobject.GetPointData().SetScalars(mapped_colors)
            end
        end
    elseif color in keys(dataset.cell_data)
        if RGB
            jl_mapped_colors = UInt8.(round.(reshape(dataset.cell_data[color], 
                (size(dataset.cell_data[color], 1), num_of_cells(dataset)))))
            mapped_colors = julia_to_vtk(jl_mapped_colors)
            mapped_colors.SetNumberOfComponents(size(dataset.cell_data[color], 1))
            vtkobject.GetCellData().SetScalars(mapped_colors)
        else
            jl_mapped_colors, cmin, cmax = get_jl_mapped_colors(vtkobject, dataset, color, component, opacity)
            if !(cmax ≈ cmin)
                mapped_colors = julia_to_vtk(jl_mapped_colors)
                vtkobject.GetCellData().SetScalars(mapped_colors)
            end
        end
    end

    if representation == :glyph
        vtkobject.GetPointData().SetActiveVectors(color)
        arrow_source = vtk.vtkArrowSource()
        glyph_filter = vtk.vtkGlyph3D()
        glyph_filter.SetSourceConnection(arrow_source.GetOutputPort())
        glyph_filter.SetScaleFactor(0.005*scale_factor)
        glyph_filter.SetColorModeToColorByScalar()
        glyph_filter.SetInputData(vtkobject)
        glyph_filter.OrientOn()
        glyph_filter.SetVectorModeToUseVector()
        glyph_filter.Update()
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputConnection(glyph_filter.GetOutputPort())
    else
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(vtkobject)
    end

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    if color == "" 
        actor.GetProperty().SetColor(1.0, 1.0, 1.0)
    end
    if isdefined(:cmin) && isdefined(:cmax) && cmin ≈ cmax
        actor.GetProperty().SetColor(0.0, 0.0, 1.0)
    end
    if representation == :wireframe
        property = vtk.vtkProperty()
        property.SetRepresentationToWireframe()
        actor.SetProperty(property)
    end

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.1, 0.2, 0.3)
    renderer.AddActor(actor)

    camera = vtk.vtkCamera()
    setup_camera(camera, dataset)
    renderer.SetActiveCamera(camera)

    renderwindow = vtk.vtkRenderWindow()
    renderwindow.SetSize(window_size...)
    renderwindow.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderwindow)

    if legend && color != "" && !RGB && !(cmin ≈ cmax)
        add_scalar_bar!(interactor, cmin, cmax, ncolors, opacity, legendtitle)
    end

    axes = vtk.vtkAxesActor() 
    widget = vtk.vtkOrientationMarkerWidget()
    widget.SetOutlineColor(0.9300,0.5700,0.1300)
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(interactor)
    widget.SetViewport(0., 0., 0.2, 0.2)
    widget.SetEnabled(1)
    widget.InteractiveOn()

    return renderwindow, interactor
end
