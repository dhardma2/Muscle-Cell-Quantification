# Muscle-Cell-Quantification
Code and macros to accompany the manuscript: “Quantitative measurement of morphometric indicators of skeletal muscle cell behaviour and quality”
Code written by David Hardman, Centre for Medical Informatics, Usher Institute, University of Edinburgh unless explicitly stated.

#Example macros used for pre-processing and/or binarising images in ImageJ 
#-> /ImageJ_Preprocessing_Macros/

#Quantification of live images of myoblast motion and fusion 
-> /Live_Images/
For stained images: CellTrackWorkFlowLiveImaging.m, for brightfield images: CellTrackWorkFlowUnStained.m
For orientation statistics on myotubes: MTAngleStats.m


#Quantification of static images 
-> /Static_Images/
For manual measurement of myonuclei and actin structures -> /Manual/
For images stained for Actin and DAPI -> Actin_DAPI_only
Un-binarised images -> StaticCellBinarizSeg.m , binarised images -> StaticCellSeg.m
For images stained for Actin, DAPI and a marker of myonuclei (Myogenin or pericentrin) -> /Actin_DAPI_myonuclei/
Option to include the average pixel area of myoblast nuclei (MCmean0) if known to improve accuracy.
If images of BTX staining to quantify acetylcholine receptors is available -> StaticCellSegBTX.m

#Additional myotube orientation metrics 
-> /Myotube_orientation_statistics/Orientation_extractor.m
