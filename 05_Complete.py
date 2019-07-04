import sys, os
sys.path.append('./pttx')
from shutil import copy2

import automatic_age_grid_seeding as aags

print 'All modules imported successfully'

##########################################################

# Set the input parameters by pointing to a yaml file

config_file = sys.argv[1]

(input_rotation_filenames, topology_features, COBterrane_file, present_day_age_grid_file,
 static_polygon_file, seedpoints_output_dir, grd_output_dir, output_gridfile_template,
 min_time, max_time, mor_time_step, gridding_time_step, ridge_sampling,
 initial_ocean_healpix_sampling, initial_ocean_mean_spreading_rate, area_threshold,
 grdspace, xmin, xmax, ymin, ymax, region, grid_masking, num_cpus) = aags.get_input_parameters(config_file)


print 'Input parameter definition completed'


if not os.path.isdir(grd_output_dir):
    os.mkdir(grd_output_dir)
if not os.path.isdir('{0}/unmasked/'.format(grd_output_dir)):
    os.mkdir('{0}/unmasked/'.format(grd_output_dir))
if not os.path.isdir('{0}/masked/'.format(grd_output_dir)):
    os.mkdir('{0}/masked/'.format(grd_output_dir))

###################################################


initial_ocean_seedpoint_filename = '%s/age_from_distance_to_mor_%0.2fMa.gmt' % (seedpoints_output_dir, max_time)
mor_seedpoint_filename = './%s/MOR_plus_one_merge_%0.2f_%0.2f.gmt' % (seedpoints_output_dir, min_time, max_time)


aags.get_initial_ocean_seeds(topology_features, input_rotation_filenames, COBterrane_file, seedpoints_output_dir,
                             max_time, initial_ocean_mean_spreading_rate, initial_ocean_healpix_sampling,
                            area_threshold, mask_sampling=grdspace)


aags.get_isochrons_from_topologies(topology_features, input_rotation_filenames, seedpoints_output_dir,
                                   max_time, min_time, mor_time_step, ridge_sampling, num_cpus=num_cpus)

aags.reconstruct_seeds(input_rotation_filenames, topology_features, seedpoints_output_dir,
                       mor_seedpoint_filename, initial_ocean_seedpoint_filename,
                       max_time, min_time, gridding_time_step, dump_to_file=True)



aags.make_grids_from_reconstructed_seeds(input_rotation_filenames, max_time, min_time, gridding_time_step,
                                         grdspace, region, grd_output_dir, output_gridfile_template,
                                         present_day_age_grid_desampling_factor=5, num_cpus=num_cpus,
                                         present_day_age_grid_file=present_day_age_grid_file,
                                         static_polygon_file=static_polygon_file,
                                         COBterrane_file=COBterrane_file)

copy2(config_file, grd_output_dir)
