import sys
#import pygplates
import os
import numpy as np
from shutil import copy2
import yaml



import automatic_age_grid_seeding as aags

print('All modules imported successfully')

##########################################################

# Set the input parameters by pointing to a yaml file

#config_file = "./config_M16a.yaml"
#config_file = "./config_M1000-NNR.yaml"
config_file = sys.argv[1]

def run_paleo_age_grids(config_file):

    (input_rotation_filenames, topology_features, COBterrane_file,
    grd_output_dir, output_gridfile_template,
    min_time, max_time, mor_time_step, gridding_time_step, ridge_sampling,
    initial_ocean_healpix_sampling, initial_ocean_mean_spreading_rate, area_threshold,
    grdspace, xmin, xmax, ymin, ymax, region, grid_masking, num_cpus, backend) = aags.get_input_parameters(config_file)



    print('Input parameter definition completed')

    subduction_collision_parameters = (5.0, 10.0)
    #continent_mask_file_pattern = None
    continent_mask_file_pattern = '%s/masks/mask_{:0.1f}Ma.nc' % grd_output_dir

    seedpoints_output_dir = '{:s}/seedpoints/'.format(grd_output_dir)

    if not os.path.isdir(grd_output_dir):
        os.mkdir(grd_output_dir)
    if not os.path.isdir('{0}/unmasked/'.format(grd_output_dir)):
        os.mkdir('{0}/unmasked/'.format(grd_output_dir))
    if not os.path.isdir('{0}/masked/'.format(grd_output_dir)):
        os.mkdir('{0}/masked/'.format(grd_output_dir))
    if not os.path.isdir('{0}/masks/'.format(grd_output_dir)):
        os.mkdir('{0}/masks/'.format(grd_output_dir))
    if not os.path.isdir('{0}/gridding_input/'.format(grd_output_dir)):
        os.mkdir('{0}/gridding_input/'.format(grd_output_dir))
    if not os.path.isdir('{0}/seedpoints/'.format(grd_output_dir)):
        os.mkdir('{0}/seedpoints/'.format(grd_output_dir))
    ###################################################


    initial_ocean_seedpoint_filename = '{:s}/seedpoints/age_from_distance_to_mor_{:0.2f}Ma.gmt'.format(grd_output_dir, max_time)
    mor_seedpoint_filename = '{:s}/seedpoints/MOR_plus_one_merge_{:0.2f}_{:0.2f}.gmt'.format(grd_output_dir, min_time, max_time)

    #'''
    aags.make_masking_grids(COBterrane_file, input_rotation_filenames, max_time, min_time, gridding_time_step,
                            grdspace, region, grd_output_dir, output_gridfile_template, num_cpus)
    #'''

    aags.get_initial_ocean_seeds(topology_features, input_rotation_filenames, COBterrane_file, seedpoints_output_dir,
                                max_time, initial_ocean_mean_spreading_rate, initial_ocean_healpix_sampling,
                                area_threshold, mask_sampling=grdspace)


    aags.get_isochrons_from_topologies(topology_features, input_rotation_filenames, seedpoints_output_dir,
                                    max_time, min_time, mor_time_step, ridge_sampling, num_cpus=num_cpus)
    #'''
    #'''
    aags.reconstruct_seeds(input_rotation_filenames, topology_features, seedpoints_output_dir,
                        mor_seedpoint_filename, initial_ocean_seedpoint_filename,
                        max_time, min_time, gridding_time_step, grd_output_dir,
                        subduction_collision_parameters=subduction_collision_parameters,
                        continent_mask_file_pattern=continent_mask_file_pattern, backend=backend)
    #'''


    aags.make_grids_from_reconstructed_seeds(input_rotation_filenames, max_time, min_time, gridding_time_step,
                                            grdspace, region, grd_output_dir, output_gridfile_template,
                                            num_cpus=num_cpus, COBterrane_file=COBterrane_file)
    #'''

    copy2(config_file, grd_output_dir)




if __name__ == "__main__":
    
    config_file = sys.argv[1]
    run_paleo_age_grids(config_file)

