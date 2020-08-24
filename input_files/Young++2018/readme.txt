Global kinematics of tectonic plates and subduction zones since the late Paleozoic Era

	

Alexander Young(1), Nicolas Flament(1), Kayla Maloney(2), Simon Williams(2), Kara Matthews(2), Sabin Zahirovic(2), Dietmar Müller(2,3)

1. The University of Wollongong, NSW 2522, Australia 

2. EarthByte Group, School of Geosciences, The University of Sydney, NSW 2006, Australia

3. Sydney Informatics Hub, The University of Sydney, NSW 2006, Australia 



Contact: ajy321@uowmail.edu.au


Supplementary Material

We provide the digital plate model files (including rotations and geometries). These files allow for the visualisation and/or manipulation of the late Paleozoic to present-day (410-0 Ma) global plate motion model presented in this study. 


#########################################
The digital plate model files are compatible with the open-source GPlates plate reconstruction software (www.gplates.org):

(1) Rotations - Global rotation model that contains the reconstruction poles that describe the motions of the continents and oceans.
- Global_250-0Ma_Young_et_al.rot (455 KB)
- Global_410-250Ma_Young_et_al.rot (154 KB)

(2) Plate polygons and boundary geometries - Topologically closed plate polygons are constructed from the intersection of ridges, transforms, subduction zones and other plate boundary geometries. These 'resolved topologies' are defined at 1 Myr intervals (410-0 Ma). The plate boundary geometries and plate polygons have been assigned plate reconstruction IDs to allow them to be reconstructed using the supplied rotation file.
- Global_Mesozoic-Cenozoic_plate_boundaries_Young_et_al.gpml (36.5 MB)
- Global_Paleozoic_plate_boundaries_Young_et_al.gpml (6.7 MB)
- TopologyBuildingBlocks_Young_et_al.gpml (2 MB) - this file is identical to Müller et al. (2016)

(3) Coastlines - Geometries of the present-day coastlines.
- Global_coastlines_Young_et_al_low_res.shp (1.2 MB  including auxiliary files, datum-WGS 1984)


(4) Static polygons (optional) - Includes ocean isochron and terrane polygon geometries.
- Global_GPlates_PresentDay_StaticPlatePolygons_Young_et_al.shp (1.4 MB inc. auxillary files, datum-WGS 1984)

(5) Continental polygons (optional) - Includes continental terrane polygon geometries and excludes oceanic lithosphere.
- PresentDay_ContinentalPolygons_Young_et_al.shp (451 KB inc. auxillary files, datum-WGS 1984)

GPlates: 
To view the model, load all files in GPlates (either drag and drop files onto the globe OR from the navigation bar at the top of the screen click File -> Open Feature Collection and select files). Both rotation files (1) and each of the three plate geometry files (2) need to be loaded for the model to work properly. It is recommended that coastlines (3) are loaded to see how the continents move, however only one coastline file is necessary (.gpml or .shp). The static polygons (4) and continental polygons (5) are optional. 

The two rotation files need to be 'connected' in order for the model to run continuously from 410 to 0 Ma. In the GPlates 'Layers' window (opened from the main navigation bar, click 'Window' -> 'Show Layers') the rotation files will be highlighted yellow, yet only one will have a yellow tick next to it to signify it is being used. Click the small black triangle to the left the ticked rotation file. Under 'Inputs' -> 'Reconstruction features' click 'Add new connection' and then select the other rotation file from the list of files that will appear. This will ensure that both rotation files are active. 

Finally, it is recommended to experiment with geometry visibility in order to make the globe less cluttered. For instance, from the navigation bar click View -> Geometry Visibility and untick 'Show Line Geometries'. Alternatively, files can be toggled on and off using the tick boxes in the Layers window. For more information about using GPlates, a set of user tutorials can be accessed from the GPlates website - http://www.gplates.org/docs.html.


#########################################
MODEL REFERENCING:
When using our model, in addition to citing this publication, please consider citing the studies of Domeier and Torsvik (2014), Matthews et al. (2016) and Müller et al. (2016) which served as the basis for this model in the late Paleozoic and Mesozoic-Cenozoic, respectively, and citing any other study that describes refinements to the plate reconstructions in your region of interest as appropriate. 

- Domeier, M., & Torsvik, T. H. (2014). Plate tectonics in the late Paleozoic. Geoscience Frontiers, 5(3), 303-350. DOI:10.1016/j.gsf.2014.01.002
- Müller, R. D., Seton, M., Zahirovic, S., Williams, S. E., Matthews, K. J., Wright, N. M., Shephard, G. E., Maloney, K., Barnett-Moore, N., Hosseinpour, M., Bower, D. J., & Cannon, J. (2016). Ocean Basin Evolution and Global-Scale Plate Reorganization Events Since Pangea Breakup. Annual Review of Earth and Planetary Sciences, 44(1). DOI:10.1146/annurev-earth-060115-012211
-Matthews, K. J., Maloney, K. T., Zahirovic, S., Williams, S. E., Seton, M., & Mueller, R. D. (2016). Global plate boundary evolution and kinematics since the late Paleozoic. Global and Planetary Change, 146, 226-250.
DOI:10.1016/j.gloplacha.2016.10.002
