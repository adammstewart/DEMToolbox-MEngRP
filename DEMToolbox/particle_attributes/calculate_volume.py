import numpy as np
import warnings

# Addition to Jack Grogan's DEMToolbox package for use in MEng RP to calculate 
# particle volumes with correction factor for multispheres based on aspect ratio.

def calculate_volume(particle_data, aspect_ratio=1.0, sphere_radius=0.0005):
    """Calculate and append the volume of particles to the particle data.

    This function calculates the volume of each particle. If an aspect ratio
    greater than 1.0 is provided, it calculates a correction factor for 
    multisphere particles based on the mathematical overlap of their subspheres.
    The corrected volume is then stored directly in the VTK point data array.

    Logic for multisphere volume correction based on multispheres being built from 
    minimum number of required subspheres to achieve specified aspect ratio 
    (e.g., aspect ratio 1.2 though 2.0 requiring 2 subspheres, aspect ratio 2.2 
    through 3.0 requiring 3 subspheres, etc.), and the total volume of the separate 
    subspheres compared to the actual known volume of the multisphere. For MEng RP
    spheres and multispheres are designed to have the same actual volume, equal to 
    that of a sphere with radius 0.0005 m.

    Parameters
    ----------
    particle_data : vtkPolyData
        The particle vtk containing the particles.
    aspect_ratio : float, optional
        The aspect ratio of the multispheres. Valid inputs are 1.0, or 
        1.2 to 3.0 in increments of 0.2. Default is 1.0 (perfect spheres).
    sphere_radius : float, optional
        The theoretical radius of a standard sphere (in meters) which shares
        the exact same actual volume as the multispheres. Default is 0.0005.

    Returns
    -------
    particle_data : vtkPolyData
        The particle vtk with the "volume" point data array added.
    correction_factor : float
        The calculated correction factor applied to the multispheres. 
        Returns 1.0 if no correction was needed.

    Raises
    ------
    KeyError
        If the 'radius' array is missing from the point data.
    ValueError
        If the provided aspect ratio is not supported by the dictionary.
    """
    if "radius" not in particle_data.point_data:
        raise KeyError("The particle data must contain a 'radius' array to calculate volume.")
    if "volume" in particle_data.point_data:
        warnings.warn("The particle data already contains a 'volume' array. It will be overwritten.", UserWarning)
        
    radii = particle_data.point_data["radius"]
    
    # 1. Base volume calculation (assuming all points are standard spheres first)
    volumes = (4.0 / 3.0) * np.pi * (radii ** 3)
    
    # Fast exit for perfect spheres
    if np.isclose(aspect_ratio, 1.0):
        particle_data.point_data["volume"] = volumes
        return particle_data, 1.0

    # 2. Stability check for integer aspect ratios
    # If the aspect ratio is a whole number (e.g., 2.0, 3.0), there is 
    # zero overlap between subspheres. The correction factor is exactly 1.0.
    if aspect_ratio % 1 == 0:
        correction_factor = 1.0
    else:
        # For non-integers, calculate the overlap correction factor
        subsphere_radii_map = {
            1.2: 0.00045872,
            1.4: 0.00043049,
            1.6: 0.00041175,
            1.8: 0.00040073,
            2.2: 0.00036446,
            2.4: 0.00035667,
            2.6: 0.00035120,
            2.8: 0.00034789
        }
    
        ar_key = round(aspect_ratio, 1)
    
        if ar_key not in subsphere_radii_map:
            raise ValueError(f"Non-integer aspect ratio {aspect_ratio} is not supported. "
                             "Please add its subsphere radius to the map.")
    
        subsphere_radius = subsphere_radii_map[ar_key]
        
        # Determine the number of subspheres using ceiling
        n_subspheres = int(np.ceil(aspect_ratio))
    
        # Calculate the mathematical correction factor
        actual_volume = (4.0 / 3.0) * np.pi * (sphere_radius ** 3)
        sum_subsphere_volume = n_subspheres * ((4.0 / 3.0) * np.pi * (subsphere_radius ** 3))
        
        correction_factor = actual_volume / sum_subsphere_volume

    # 3. Apply the correction factor specifically to multisphere elements
    # Identify multispheres as any particle whose radius doesn't match the sphere_radius
    multisphere_mask = ~np.isclose(radii, sphere_radius)

    if not np.any(multisphere_mask) and aspect_ratio != 1.0:
        warnings.warn(f"Aspect ratio {aspect_ratio} was specified, but no multispheres "
                      f"(radii differing from {sphere_radius}) were found in the data. "
                      "Correction factor calculated but not applied.", UserWarning)
        
    volumes[multisphere_mask] *= correction_factor

    # 4. Assign the final corrected volumes back to the VTK dataset
    particle_data.point_data["volume"] = volumes

    return particle_data, correction_factor