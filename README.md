# MSThesis
Intensity Augmented ICP (IAICP)

Automatic registration of Laser Scanner data is of great interest in Geoinformatics research. Conventional LiDAR data registration techniques require much human intervention in the field and at data processing stage. This dissertation presents a novel approach for Laser Scanner data registration utilizing radiometric data collected by the scanner along with the geometric coordinates for estimating the transformation parameters between the unregistered point clouds. Before using radiometric data in registration process the thesis proposes to normalize these.  This thesis, therefore, also develops a new technique for the normalization of LiDAR intensity information so that the intensity data can be effectively incorporated along with the geometric attributes to perform registration. The proposed registration algorithm is a variant of the conventional Iterative Closest Point (ICP) algorithm which is the classical solution for all types of 3D registration problems. The algorithm presented in this thesis works in two stages with Intensity Augmented ICP (IAICP) for coarse registration stage and conventional geometric ICP at the fine registration stage. The proposed approach is successfully applied to a set of benchmarked data resulting in an accurate estimation of the transformation parameters. A comparison of the conventional and intensity augmented registration approaches is also presented.  The results obtained for intensity normalization indicate that the proposed scheme results in bringing two different scans of same object closer in terms of their radiometric characteristics.  Further results indicate the supremacy of IAICP over the ICP, as the latter is found to fail in geometrically confusing cases while intensity augmented ICP works in these cases.  

Cite:

An evaluation of intensity augmented ICP for terrestrial LiDAR data registration 
Journal of Geomatics, Vol. XI No.2, pp 139-148 October, 2017 Authors: Bharat Lohani, Sandeep Sasidharan

Intensity Augmented ICP for Registration of Laser Scanner Point Clouds
Proceedings of XXXII INCA International Congress on Cartography for Sustainable Earth Resource Management, Indian Cartographer, Vol. XXXII, pp 30-34 September 17, 2013 Authors: Bharat Lohani, Sandeep Sasidharan

A Normalization scheme for Terrestrial LiDAR Intensity Data by Range and Incidence Angle 
International Journal of Emerging Technology and Advanced Engineering, Vol. 6, Issue. 5, pp 322-328 May, 2016 Author: Sandeep Sasidharan

Intensity Augmented ICP for 3D Laser Scanner Point Cloud Registration
Master Thesis, Sandeep Sasidharan, IIT Kanpur 2012

#Code reference
https://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point?focused=3775639&tab=function
