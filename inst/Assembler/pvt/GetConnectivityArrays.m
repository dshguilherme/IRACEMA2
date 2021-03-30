function [global_basis_index, element_local_mapping, element_ranges] = ...
    GetConnectivityArrays(GeometryObj)
    switch GeometryObj.rank
        case 1
            [global_basis_index, element_local_mapping, element_ranges] ...
                = CurveConnectivity(GeometryObj);
        case 2
            [global_basis_index, element_local_mapping, element_ranges] ...
                = SurfaceConnectivity(GeometryObj);
        case 3
            [global_basis_index, element_local_mapping, element_ranges] ...
                = VolumeConnectivity(GeometryObj);
    end
end