# agegrid-dev
Development repo for auto-age-gridding workflow

Repository contains a workflow to generate age grids from a topological plate reconstruction

The entire workflow can be launched using the python file '05_Complete.py'

The inputs are contained in a config file <filename>.yaml, where you enter all the necessary input parameters (input file names and directories, values that define the age range, grid resolution, etc...)

## Requirements

The following python libraries are required:
- pygplates
- numpy
- pandas
- xarray
- joblib
- PlateTectonicsTools

The code also calls GMT functions (blockmedian, sphinterpolate, grdmask, gmtselect, grdtrack) assuming a GMT version 5 syntax.

TODO - reorganise dependencies

# Workflow Structure

### get_input_parameters
*Retrieves the input parameters that control the seafloor reconstruction process from yaml file*

An example input file listing all parameters looks like this (in this case using the GPlates 2.1 sample data).

    MODELDIR: './input_files/'

        input_rotation_filenames:
            - 'Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot'

        topology_features:
            - 'DynamicPolygons/Matthews_etal_GPC_2016_MesozoicCenozoic_PlateTopologies.gpml'
            - 'DynamicPolygons/Matthews_etal_GPC_2016_Paleozoic_PlateTopologies.gpml'
            - 'DynamicPolygons/Matthews_etal_GPC_2016_TopologyBuildingBlocks.gpml'

        COBterrane_file: 'StaticPolygons/Global_EarthByte_GeeK07_COB_Terranes_Matthews_etal.gpml'

        present_day_age_grid_file: None
        static_polygon_file: None

    # name of output directories
    OutputFiles:
        seedpoints_output_dir: "./seed_point_files/"
        grd_output_dir: "./tmp/"
        #grd_output_dir: './velocity_scaling/test/'
        output_gridfile_template: 'M16_seafloor_age_'

    TimeParameters:
    # --- set parameters
        min_time: 0.
        max_time: 410.
        mor_time_step: 1.
        gridding_time_step: 1.

    # Distance in arc-degrees along ridges at which to tessellate lines to
    # create seed points (note: does not desample the points where the 
    # original line geometries
    # are already closer than the specified distance
    SpatialParameters:
        ridge_sampling: 1.0
        initial_ocean_healpix_sampling: 32
        initial_ocean_mean_spreading_rate: 75.

        # used to remove small polygons in the masking. Units are area on unit sphere
        area_threshold: 0.0001

        seed_age_offset: 0.0

        # Control the gridding extent and resolution parameters passed to GMT
        grdspace: 0.25
        xmin: -180.
        xmax: 180.
        ymin: -90.
        ymax: 90.

        grid_masking: True

    num_cpus: 12


### get_initial_ocean_seeds
*Generates a set of points at the oldest time selected for the reconstruction sequence.*

The points populate the ocean basins (points are masked out within the areas defined by the ‘COB-Terranes’ polygons) and are assigned ages assuming a uniform average spreading rate combined with the distance of each point to the nearest feature of type ‘MidOceanRidge’ feature in the resolved plate boundary of the the plate containing the point. Note that the approach taken is fairly simplistic and does not (yet) take account of where continents are located that may suggest the nearest ridge to a given point was not ridge where that section of seafloor was created.

The resulting points are saved into a single file.
get_isochrons_from_topologies
Generates points along the features of type ‘MidOceanRidge’ extracted from a time-sequence of resolved topological plate boundary configurations

The code here generates resolved topologies at a defined interval (defined in Myr). This number can be defined as any value including fractions of Myr, but keep in mind that the topological reconstructions are most likely to be robust (without topological error) at integer values of time. 

For each resolved topology, the MidOceanRidge features are isolated (other types of plate boundary are ignored). These geometries include both ridge and transform segments. The ridge segments are isolated using the function separate_ridge_transform_segments contained within the package ‘PlateTectonicsTools’ (https://github.com/EarthByte/PlateTectonicTools). Seed points are only generated along these ridge segments.

The density of the points generated along each ridge segment is defined by the parameter ‘ridge_sampling’. Note that this sampling is used as a ‘target’ sampling, but the original nodes of the lines are always preserved, hence the actual point sampling will always be more dense that the specified value (and will vary from one segment to the next). See the pygplates documentation on polyline tesselation for more detail.

To allow the points to be assigned a plate id for the ‘reconstruct by topologies’ stage, the points lying along the ridge geometries are duplicated. These two points are shifted a small amount (set to 0.01 degrees) away from the ridge (one point in each direction) based on the spreading direction (itself derived from the stage pole describing the opening between the relevant plates). 

The resulting points are saved to a file, one file per time increment.

### reconstruct_seeds
*Combines the points generated from both the ‘get_initial_ocean_seeds’ and ‘get_isochrons_from_topologies’ steps, and reconstructs them using the ‘reconstruct_by_topologies’ method.*

The concept of reconstructing feature geometries using topologies was introduced by Müller et al. (2018). The code used to generate seafloor age distributions uses this concept. Some parameters must be defined to handle cases where the plate id of one or many points changes between two successive time steps, so that cases where points have been subducted (and should therefore be deactivated) can be distinguished from points that should persist following a plates splitting or merging within an ocean basin.
 
The following text is a detailed description of the logic applied to each point at each time step, if the point is found to be in a different plate than the previous time step (see also the code comments within the file ‘reconstruct_by_topologies.py’: 

If transitioning from a rigid plate to another rigid plate with a different plate ID then calculate the difference in velocities and continue testing as follows (otherwise, if there's no transition, then the point is still active)...

If the velocity difference is below a threshold then we assume the previous plate was split, or two plates joined. In this case the point has not subducted (forward in time) or been consumed by a mid-ocean (backward in time) and hence is still active.

If the velocity difference is large enough then we see if the distance of the *previous* position to the polygon boundary (of rigid plate containing it) exceeds a threshold. If the distance exceeds the threshold then the point is far enough away from the boundary that it cannot be subducted or consumed by it and hence the point is still active. However if the point is close enough then we assume the point was subducted/consumed (remember that the point switched plate IDs). Also note that the threshold distance increases according to the velocity difference to account for fast moving points (that would otherwise tunnel through the boundary and accrete onto the other plate). The reason for testing the distance from the *previous* point, and not from the *current* point, is:

   (i)  A topological boundary may *appear* near the current point (such as a plate split at the current time) and we don't want that split to consume the current point regardless of the velocity difference. It won't get consumed because the *previous* point was not near a boundary (because before split happened). If the velocity difference is large enough then it might cause the current point to transition to the adjacent split plate in the *next* time step (and that's when it should get consumed, not in the current time step). An example of this is a mid-ocean ridge suddenly appearing (going forward in time).

   (ii) A topological boundary may *disappear* near the current point (such as a plate merge at the current time) and we want that merge to consume the current point if the velocity difference is large enough. In this case the *previous* point is near a boundary (because before plate merged) and hence can be consumed (provided velocity difference is large enough). And since the boundary existed in the previous time step, it will affect position of the current point (and whether it gets consumed or not). An example of this is a mid-ocean ridge suddenly disappearing (going backward in time).

[note that items (i) and (ii) above apply both going forward and backward in time.]

The parameters that control the collision detection are therefore:
threshold velocity delta in kms/my:  default value is 7 kms/my
threshold distance to boundary per My in kms/my:  default value is 10 kms/yr

These parameters can be defined individually for different types of plate boundary using the ‘feature_specific_collision_parameters’ option. For example, we may want to use more aggressive rejection parameters when detecting collisions near to features that are labelled as type ‘SubductionZone’.

The results of this step are that the reconstructed points at any given snapshot are saved together into a single file (ascii format) for that snapshot with coordinates and seafloor age at the time, suitable to be interpolated onto a regular grid.

### make_grids_from_reconstructed_seeds
*Interpolate reconstructed points to regular grids stored in netcdf format*

The ascii files generated by the reconstruct_seeds function are transformed from their irregular sampling locations to a regular grid in geographich coordinates using two GMT functions. First, the points are passed through the function blockmedian to create a set of points with a spatial sampling similar to the desired grid spacing. These points are then gridded using piecewise linear interpolation in spherical coordinates as implemented in the GMT function sphinterpolate. 
