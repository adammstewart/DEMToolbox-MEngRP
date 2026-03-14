import numpy as np
import warnings

from ..classes.particle_samples import ParticleSamples
from ..classes.particle_attribute import ParticleAttribute

# Added to Jack Grogan's DEMToolbox for MEng RP for sampling by shape (sphere vs multisphere)
 

def sample_by_shape(particle_data, 
                    sphere_radius=0.0005, 
                    append_column="shape_sample",
                    particle_id_column="id"):
    """Sample the particles based on whether they are spheres or multispheres.

    This function classifies particles into two distinct samples: 
    0 (Spheres) and 1 (Multisphere subspheres). The classification is determined by 
    comparing the particle's radius to the specified sphere_radius. The 
    resulting particle data will have a new column added with the sample 
    class for each particle, and a ParticleSamples object will be returned.

    Sampling logic classifies using particle radius, and as such spheres and 
    multispheres must have different radii. For MEng RP, all spheres had a radius
    of 0.0005 m, while subspheres have smaller radii depending on required aspect 
    ratio of the multisphere. 

    Parameters
    ----------
    particle_data : vtkPolyData
        The particle vtk containing the particles to be sampled.
    sphere_radius : float, optional
        The theoretical radius of a standard sphere. Particles with a radius
        not equal to this value are classified as multispheres. Default is 0.0005.
    append_column : str, optional
        The name of the samples column to append to the particle data.
        Default is "shape_sample".
    particle_id_column : str, optional
        The name of the particle id column in the particle data, by
        default "id".

    Returns
    -------
    particle_data : vtkPolyData
        The particle vtk with the attribute column added.
    samples : ParticleSamples
        A ParticleSamples object containing the samples column name,
        the particle ids and their corresponding sample ids stored in a
        ParticleAttribute object, a list of sample elements (0 and 1), 
        a list of occupied sample elements, the number of particles in
        each sample, and the total sampled and unsampled particles.
        
    Raises
    ------
    KeyError
        If the 'radius' or particle id column is missing from the point data.
    UserWarning
        If the particle data has no points, a warning is issued and
        an empty ParticleSamples object is returned.
    """
    
    # 1. Handle empty datasets
    if particle_data.n_points == 0:
        warnings.warn("Cannot sample empty particles file", UserWarning)
        sample_attribute = ParticleAttribute(particle_id_column, 
                                             append_column,
                                             np.empty((0, 2)))
        samples = ParticleSamples(
            append_column, sample_attribute, [], [], [], 0, 0)
        return particle_data, samples
        
    # 2. Check for required data arrays
    if "radius" not in particle_data.point_data:
        raise KeyError("The particle data must contain a 'radius' array to sample by shape.")
    if particle_id_column not in particle_data.point_data:
        raise KeyError(f"The particle data must contain a '{particle_id_column}' array.")

    # 3. Classify the particles
    radii = particle_data.point_data["radius"]
    
    # ~np.isclose returns True for multispheres, False for standard spheres
    multisphere_mask = ~np.isclose(radii, sphere_radius)
    
    # Convert the boolean mask to integers: 0 (Sphere) and 1 (Multisphere)
    samples_column = multisphere_mask.astype(int)

    # 4. Gather sample statistics
    # Our possible cells are always 0 and 1
    cells = np.array([0, 1], dtype=int) 
    
    occupied_cells, counts = np.unique(samples_column, return_counts=True)
    
    cell_particles = np.zeros(2, dtype=int)
    np.add.at(cell_particles, occupied_cells, counts)
    
    n_sampled_particles = particle_data.n_points
    n_unsampled_particles = 0 # In a binary shape check, all particles get sampled

    # 5. Add the samples column to the particle data
    particle_data.point_data[append_column] = samples_column

    # 6. Create the attribute and sample objects
    sample_data = np.array([particle_data.point_data[particle_id_column], samples_column]).T
    sample_attribute = ParticleAttribute(particle_id_column, 
                                         append_column,
                                         sample_data)

    # Note: Removed vector_1_centers and vector_1_bounds as they aren't relevant for categorical sampling
    samples = ParticleSamples(append_column,
                              sample_attribute,
                              cells,
                              occupied_cells,
                              cell_particles,
                              n_sampled_particles,
                              n_unsampled_particles)

    return particle_data, samples