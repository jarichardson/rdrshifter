#include "structs.h"

int main(int argc, char *argv[]) {
	
	/*VARIABLES******************************************************************/
	/*Files*/
	char     *configFilename = argv[1]; /*configuration file path               */
	char     **Filenames;               /*A list of file paths for model output */
	
	point *DEMvals;
	point *RDRvals;
	point *shftRDRvals;
	point *cropDEMvals;
	point  origin;
	point  translate;
	point  rotate;
	double scale;
	
	point dist;
	double hypotenuse;
	double curangle;
	point *RDRinterp;
	point RDRUL;
	point RDRLR;
	point DEMUL;
	point DEMLR;	
	point expandedDEMUL;
	point expandedDEMLR;
	unsigned RDRinboundsct=0;
	unsigned DEMinboundsct=0;
	
	double NNdist[4];
	double NNelev[4];
	double NNwt[4];
	double curdistance;
	
	double DEMnoDATA;
	double RDRnoDATA;
	unsigned demmax;
	unsigned rdrmax;
	
	double *DEMmetadata;
	double *RDRmetadata;
	
	FILE *DEMfile;
	FILE *RDRfile;
	
	int ret, i, j;
	
	/*TIME **********************************************************************/
	startTime = time(NULL); /*Define Start Time*/
	
	/*WELCOME USER AND CHECK FOR CORRECT USAGE***********************************/
	printf("\n\n               This will shift ur data.\n\n");
	
	/*User must supply the name of the executable and a configuration file*/
	if(argc<2) {
		printf("Usage: %s config-filename\n",argv[0]);
		return(-1);
	}
	
	printf("Reading Configuration File...\n");
	
	/*MODULE: INITIALIZE*********************************************************/
	/*        Assigns several empty variables based on a user-defined 
	            configuration file.                                             */
	/*        File Name List output in this order:
	            [0] - DEM
	            [1] - Residual Flow Thickness
	            [2] - Elevation Uncertainty
	            [3] - Output file: ASCII X,Y,Thickness
	            [4] - Output file: Hit Map
	            [5] - Output file: Raster Thickness
	            [6] - Output file: Raster Elevation
	            [7] - Output file: Raster Elevation + Flow Thickness            */
	
	ret = INITIALIZE(configFilename,        /*chr Configuration File Name       */
	                 &Filenames,            /*chr File Name List                */
	                 &origin.easting,       /*point  origin location            */
	                 &origin.northing,       /*point  origin location            */
	                 &origin.elevation,       /*point  origin location            */
	                 &translate.easting,       /*point  origin location            */
	                 &translate.northing,       /*point  origin location            */
	                 &translate.elevation,       /*point  origin location            */
	                 &rotate.easting,       /*point  origin location            */
	                 &rotate.northing,       /*point  origin location            */
	                 &rotate.elevation,       /*point  origin location            */
	                 &scale,                /*double scale                      */
	                 &DEMnoDATA,            /*DEM raster No Data                */
	                 &RDRnoDATA             /*RDR raster No Data                */
	                );

	/*Check for Error flag (INITIALIZE returns <0 value)*/
	if(ret<0){
		printf("\nERROR [MAIN]: Error flag returned from [INITIALIZE].\n");
		printf("Exiting.\n");
		return(-1);
	}
	
	printf("Loading Digital Elevation Model...");
	/*MODULE: DEM_LOADER*********************************************************/
	/*        Loads Raster into a list and reports range values*/
	
	/*Load topography from "Expected" DEM*/
	DEMmetadata = DEM_LOADER(Filenames[0], /*char         DEM file name */
	                         &DEMvals      /*point        DEM x,y,z list*/
	                        );
	/*Check for Error flag (DEM_LOADER returns a null metadata list)*/
	if(DEMmetadata==NULL){
		printf("\nError [MAIN]: Error flag returned from DEM_LOADER [DEM].\n");
		printf("Exiting.\n");
		return(-1);
	}
	
	demmax = DEMmetadata[2]*DEMmetadata[4];
	
	/*DEMGeoTransform[0] lower left x
	  DEMGeoTransform[1] w-e pixel resolution
	  DEMGeoTransform[2] number of cols, assigned manually in this module 
	  DEMGeoTransform[3] lower left y
	  DEMGeoTransform[4] number of lines, assigned manually in this module
	  DEMGeoTransform[5] n-s pixel resolution (negative value) */
	
	DEMUL.easting = DEMmetadata[0];
	DEMUL.northing =DEMmetadata[3] + (DEMmetadata[4]*DEMmetadata[5]);
	DEMLR.easting = DEMmetadata[0] + (DEMmetadata[2]*DEMmetadata[1]);
	DEMLR.northing = DEMmetadata[3];
	expandedDEMUL.easting = DEMUL.easting - (2*DEMmetadata[1]);
	expandedDEMUL.northing = DEMUL.northing + (2*DEMmetadata[5]);
	expandedDEMLR.easting = DEMLR.easting + (2*DEMmetadata[1]);
	expandedDEMLR.northing = DEMLR.northing - (2*DEMmetadata[5]);

	printf("Loading Radar Elevation Model...");
	/*Load topography from "observed" Radar DEM*/
	RDRmetadata = DEM_LOADER(Filenames[1], /*char         RDR file name */
	                         &RDRvals      /*point        RDR x,y,z list*/
	                        );
	/*Check for Error flag (DEM_LOADER returns a null metadata list)*/
	if(RDRmetadata==NULL){
		printf("\nError [MAIN]: Error flag returned from DEM_LOADER [RDR].\n");
		printf("Exiting.\n");
		return(-1);
	}
	
	rdrmax = RDRmetadata[2]*RDRmetadata[4];
	
	RDRUL.easting = RDRmetadata[0];
	RDRUL.northing =RDRmetadata[3] + (RDRmetadata[4]*RDRmetadata[5]);
	RDRLR.easting = RDRmetadata[0] + (RDRmetadata[2]*RDRmetadata[1]);
	RDRLR.northing = RDRmetadata[3];
	
	
	/*Print range of DEM and RADAR coordinates, pre-shift.*/
	printf("DEM (Lab)    Range is:   %0.3f\n",   DEMUL.northing);
	printf("                  %0.3f    %0.3f\n", DEMUL.easting, DEMLR.easting);
	printf("                         %0.3f\n\n",   DEMLR.northing);
	
	printf("Radar (Body) Range is:   %0.3f\n",   RDRUL.northing);
	printf("                  %0.3f    %0.3f\n", RDRUL.easting, RDRLR.easting);
	printf("                         %0.3f\n\n",   RDRLR.northing);
	
	
	
	
	/*****DOF-SHIFT******/
	printf("\n\nAPPLYING ROTATIONS AND TRANSFORMATIONS TO THE LAB DATA SET\n\n");
	printf("Origin Location:\n");
	printf("         easting:   %0.3f\n" , rotate.easting);
	printf("         northing:  %0.3f\n" , rotate.northing);
	printf("         elevation: %0.3f\n" , rotate.elevation);
	
	
	/*Change Rotation from degrees to radians*/
	printf("\nRotating (CCW)...\n");
	printf("         about the Z-axis,  %0.3f° " , rotate.elevation);
	rotate.elevation *= M_PI / 180.0;
	printf("(%0.3f radians)\n", rotate.elevation);
	
	printf("         about the NS-axis, %0.3f° " , rotate.northing);
	rotate.northing   *= M_PI / 180.0;
	printf("(%0.3f radians)\n", rotate.northing);
	
	printf("         about the EW-axis, %0.3f° " , rotate.easting);
	rotate.easting  *= M_PI / 180.0;
	printf("(%0.3f radians)\n", rotate.easting);
	
	
	printf("\nTranslating...\n");
	printf("         %0.3f m to the East\n",translate.easting);
	printf("         %0.3f m to the North\n",translate.northing);
	printf("         %0.3f m        Up\n",translate.elevation);
	printf("\nScaling...\n");
	printf("         %0.3f%% of original size\n\n",(scale*100.0));
	
	for(i=0;i<rdrmax;i++){
		if (RDRvals[i].elevation != RDRnoDATA) { /*only shift if this location has a value*/
		
			/*Remove origin from locations*/
			RDRvals[i].easting   -= origin.easting;
			RDRvals[i].northing  -= origin.northing;
			RDRvals[i].elevation -= origin.elevation;
			if(i==0) {
				printf("pre-rot:  %0.3fE\t%0.3fN\t%0.3fElev\n",RDRvals[i].easting,RDRvals[i].northing,RDRvals[i].elevation);
			}
		
			/*ROTATE RDR locations**********************/
			/*ABOUT Z AXIS*/
			hypotenuse = pow((pow(RDRvals[i].easting,2.0) + pow(RDRvals[i].northing,2.0)),0.5);
			if(RDRvals[i].easting    == 0.0) {
				curangle = M_PI/2.0;
				if(RDRvals[i].northing  < 0.0) curangle *= -1.0;
			}
			else {
				curangle = atan(RDRvals[i].northing / RDRvals[i].easting);
				if(RDRvals[i].easting < 0.0) curangle += M_PI;
			}
		
			/*apply rotation*/
			RDRvals[i].easting  = cos(rotate.elevation + curangle) * hypotenuse;
			RDRvals[i].northing = sin(rotate.elevation + curangle) * hypotenuse;
		
			if(i==0) {
				printf("postrotZ: %0.3fE\t%0.3fN\t%0.3fElev\n",RDRvals[i].easting,RDRvals[i].northing,RDRvals[i].elevation);
			}
		
			/*ABOUT Y AXIS*/		
			hypotenuse = pow((pow(RDRvals[i].easting,2.0) + pow(RDRvals[i].elevation,2.0)),0.5);
			if(RDRvals[i].elevation    == 0.0) {
				curangle = M_PI/2.0;
				if(RDRvals[i].easting < 0.0) curangle *= -1.0;
			}
			else {
				curangle = atan(RDRvals[i].easting / RDRvals[i].elevation);
				if(RDRvals[i].elevation < 0.0) curangle += M_PI;
			}
		
			/*apply rotation*/
			RDRvals[i].elevation = cos(rotate.northing + curangle) * hypotenuse;
			RDRvals[i].easting   = sin(rotate.northing + curangle) * hypotenuse;
		
			if(i==0) {
				printf("postrotY: %0.3fE\t%0.3fN\t%0.3fElev\n",RDRvals[i].easting,RDRvals[i].northing,RDRvals[i].elevation);
			}
		
			/*ABOUT X AXIS*/
			hypotenuse = pow((pow(RDRvals[i].elevation,2.0) + pow(RDRvals[i].northing,2.0)),0.5);
			if(RDRvals[i].northing    == 0.0) {
				curangle = M_PI/2.0;
				if(RDRvals[i].elevation < 0.0) curangle *= -1.0;
			}
			else {
				curangle = atan(RDRvals[i].elevation / RDRvals[i].northing);
				if(RDRvals[i].northing < 0.0) curangle += M_PI;
			}
		
			/*apply rotation*/
			RDRvals[i].northing  =  cos(rotate.easting + curangle) * hypotenuse;
			RDRvals[i].elevation =  sin(rotate.easting + curangle) * hypotenuse;
		
			if(i==0) {
				printf("postrotX: %0.3fE\t%0.3fN\t%0.3fElev\n",RDRvals[i].easting,RDRvals[i].northing,RDRvals[i].elevation);
			}
		
		
			/*find distance from origin (0,0)*/
			dist.easting   = RDRvals[i].easting;
			dist.northing  = RDRvals[i].northing;
			dist.elevation = RDRvals[i].elevation;
		
			/*scale RDR locations (scale before translation)*/
			dist.easting   *= scale - 1.0; /*This is the scale shift*/
			dist.northing  *= scale - 1.0;
			dist.elevation *= scale - 1.0;
		
			/*Add translation, scale, and readd origin at same time*/
			RDRvals[i].easting   += translate.easting   + dist.easting   + origin.easting;
			RDRvals[i].northing  += translate.northing  + dist.northing  + origin.northing;
			RDRvals[i].elevation += translate.elevation + dist.elevation + origin.elevation;
		
			if(i==0) {
				printf("post-trn: %0.3fE\t%0.3fN\t%0.3fElev\n",RDRvals[i].easting,RDRvals[i].northing,RDRvals[i].elevation);
			}
		}
	}
	
	j=0;
	for(i=0;i<rdrmax;i++) {
		if (RDRvals[i].elevation != RDRnoDATA) { /*only count if this location has a value*/
			/*Find New Range*/
			if((j++)==0) {
				RDRUL.easting  = RDRvals[i].easting;
				RDRUL.northing = RDRvals[i].northing;
				RDRLR.easting  = RDRvals[i].easting;
				RDRLR.northing = RDRvals[i].northing;
			}
			else {
				if (RDRUL.easting > RDRvals[i].easting) RDRUL.easting = RDRvals[i].easting;
				else if (RDRLR.easting < RDRvals[i].easting) RDRLR.easting = RDRvals[i].easting;
				if (RDRUL.northing < RDRvals[i].northing) RDRUL.northing = RDRvals[i].northing;
				else if (RDRLR.northing > RDRvals[i].northing) RDRLR.northing = RDRvals[i].northing;
			}
		
			/*Count number of RADAR points now in or near DEM*/
			if ((RDRvals[i].easting >= expandedDEMUL.easting) && (RDRvals[i].easting <= expandedDEMLR.easting)){
				if ((RDRvals[i].northing <= expandedDEMUL.northing) && (RDRvals[i].northing >= expandedDEMLR.northing)) {
					++RDRinboundsct; /*Increment count of points within the expanded DEM boundaries*/
				}
			}
		}
	}
	
	/*Print new range of RADAR coordinates, post shift.*/
	printf("\nNew Radar Range is:   %0.3f\n",   RDRUL.northing);
	printf("               %0.3f    %0.3f\n", RDRUL.easting, RDRLR.easting);
	printf("                      %0.3f\n",   RDRLR.northing);
	printf("%u radar locations with data are within 2 pixels of the DEM\n\n", RDRinboundsct);
	
	if(RDRinboundsct==0){
		printf("NO shifted radar points lie within the DEM region!! unable to fit.\nexiting.\n");
		return(-1);
	}
	
	/*CROP ALL THE DATA*************************************************************************/
	printf("Cropping the BODY data to the LAB data extent\n");
	
	/*move Radar values to new, potentially smaller list, for speed*/
	if((shftRDRvals = malloc (sizeof (point) * (RDRinboundsct)))==NULL) {
		printf("ERROR [MAIN]: Out of Memory creating Shifted RADAR Values List!\n");
		return(-1);
	}
	
	j=0;
	for(i=0;i<rdrmax;i++){
		if((RDRvals[i].easting >= expandedDEMUL.easting) && (RDRvals[i].easting <= expandedDEMLR.easting)){
			if((RDRvals[i].northing <= expandedDEMUL.northing) && (RDRvals[i].northing >= expandedDEMLR.northing)) {
				if (RDRvals[i].elevation != RDRnoDATA) { /*only copy if there's data*/
					shftRDRvals[j] = RDRvals[i];
					j++;
				}
			}
		}
	}
	
	free(RDRvals); /*free memory of old array*/
	printf("                      copied %u values to cropped BODY array\n",(j+1));
	
	printf("Cropping the LAB data to the BODY data extent\n");
	
	/*Now find new range of DEM and scrap values outside of the Radar range*/
	if(DEMUL.easting < RDRUL.easting) DEMUL.easting = RDRUL.easting;
	if(DEMLR.easting > RDRLR.easting) DEMLR.easting = RDRLR.easting;
	if(DEMUL.northing > RDRUL.northing) DEMUL.northing = RDRUL.northing;
	if(DEMLR.northing < RDRLR.northing) DEMLR.northing = RDRLR.northing;
	
	DEMinboundsct=0;
	for(i=0;i<demmax;i++){
		if (DEMvals[i].easting >= DEMUL.easting && DEMvals[i].easting <= DEMLR.easting){
			if (DEMvals[i].northing <= DEMUL.northing && DEMvals[i].northing >= DEMLR.northing) {
				DEMinboundsct++;
			}
		}
	}
	
		/*DEMGeoTransform[0] lower left x
	  DEMGeoTransform[1] w-e pixel resolution
	  DEMGeoTransform[2] number of cols, assigned manually in this module 
	  DEMGeoTransform[3] lower left y
	  DEMGeoTransform[4] number of lines, assigned manually in this module
	  DEMGeoTransform[5] n-s pixel resolution (negative value) */

	/*move DEM values to new, potentially smaller list, for speed*/
	if((cropDEMvals = malloc (sizeof (point) * (DEMinboundsct)))==NULL) {
		printf("ERROR [MAIN]: Out of Memory creating Cropped DEM Values List!\n");
		return(-1);
	}
	
	j=0;
	for(i=0;i<demmax;i++){
		if (DEMvals[i].easting >= DEMUL.easting && DEMvals[i].easting <= DEMLR.easting){
			if (DEMvals[i].northing <= DEMUL.northing && DEMvals[i].northing >= DEMLR.northing) {
				cropDEMvals[j] = DEMvals[i];
				j++;
			}
		}
	}
	
	free(DEMvals); /*free memory of old array*/
	printf("                      copied %u values to cropped LAB array\n",(j+1));
	
	
	if((RDRinterp = malloc (sizeof (point) * (DEMinboundsct)))==NULL) {
		printf("ERROR [MAIN]: Out of Memory creating Interpolated RADAR Array!\n");
		return(-1);
	}
	
	/*Interpolate RDR onto SRTM***************************************************/
	printf("\n\nInterpolating shifted points at DEM locations\n");
	
	/*initialize NNelevs to avoid make warning.*/
	NNelev[0] = NNelev[1] = NNelev[2] = NNelev[3] = 0.0;
	
	for(i=0;i<DEMinboundsct;i++){
		/*printf("\r%u/%u",i,DEMinboundsct);*/
		
		RDRinterp[i].easting  = cropDEMvals[i].easting;
		RDRinterp[i].northing = cropDEMvals[i].northing;
		RDRinterp[i].elevation = -9999; /*starting elevation in case no value can be assigned*/
		
		/*if DEM is not NoDATA*/
		if (cropDEMvals[i].elevation != DEMnoDATA) {
		
			/*use near neighbor to find elevation
				loop through all, if in a quadrant, check if it's the closest point.
				after loop, closest 4 points get averaged, weighted by distance.
				must have all 4 points. All that is needed is to preserve a distance
				and an elevation for 4 quadrants. double[4] distance, double[4] elev*/
		
			for(j=0;j<=3;j++) {
				/*reset distance*/
				NNdist[j] = DBL_MAX;
			}
		
		
			for(j=0;j<RDRinboundsct;j++){
				if (abs(shftRDRvals[j].easting-RDRinterp[i].easting) < (2.0*DEMmetadata[1])) {
				if (abs(shftRDRvals[j].northing-RDRinterp[i].northing) < (2.0*DEMmetadata[1])) {
					curdistance = pow(pow((shftRDRvals[j].easting  - RDRinterp[i].easting ),2.0) +
					              pow((shftRDRvals[j].northing - RDRinterp[i].northing),2.0),0.5);
					
					/*First quadrant*/
					if((shftRDRvals[j].easting  >= RDRinterp[i].easting) && 
						 (shftRDRvals[j].northing >= RDRinterp[i].northing)) {
						/*is distance smallest?*/
						if (curdistance < NNdist[0]) {
							NNdist[0]  = curdistance;
							NNelev[0]  = shftRDRvals[j].elevation;
						}
					}
					/*Second quadrant*/
					else if((shftRDRvals[j].easting   < RDRinterp[i].easting) && 
						 (shftRDRvals[j].northing >= RDRinterp[i].northing)) {
						if (curdistance < NNdist[1]) {
							NNdist[1]  = curdistance;
							NNelev[1] = shftRDRvals[j].elevation;
						}
					}
					/*Third quadrant*/
					else if((shftRDRvals[j].easting  < RDRinterp[i].easting) && 
						 (shftRDRvals[j].northing < RDRinterp[i].northing)) {
						if (curdistance < NNdist[2]) {
							NNdist[2]  = curdistance;
							NNelev[2] = shftRDRvals[j].elevation;
						}
					}
					/*Fourth quadrant*/
					else {
						if (curdistance < NNdist[3]) {
							NNdist[3]  = curdistance;
							NNelev[3] = shftRDRvals[j].elevation;
						}
					}
				}} /*If the point is within 2 pixels of the point*/
			} /*For all radar values*/
		
			/*if all distances are assigned*/
			if (((NNdist[0]<DBL_MAX)&&(NNdist[1]<DBL_MAX))&&
				  ((NNdist[2]<DBL_MAX)&&(NNdist[3]<DBL_MAX))) {
				/*interpolate elevation*/
				if (NNdist[0]==0) RDRinterp[i].elevation = NNelev[0];
				else if (NNdist[1]==0) RDRinterp[i].elevation = NNelev[1];
				else if (NNdist[2]==0) RDRinterp[i].elevation = NNelev[2];
				else if (NNdist[3]==0) RDRinterp[i].elevation = NNelev[3];
				else {
					NNwt[0] = pow((1.0/NNdist[0]),2.0);
					NNwt[1] = pow((1.0/NNdist[1]),2.0);
					NNwt[2] = pow((1.0/NNdist[2]),2.0);
					NNwt[3] = pow((1.0/NNdist[3]),2.0);
					RDRinterp[i].elevation = 	(NNwt[0]*NNelev[0] + NNwt[1]*NNelev[1] + 
					                           NNwt[2]*NNelev[2] + NNwt[3]*NNelev[3]) /
					                          (NNwt[0] + NNwt[1] + NNwt[2] + NNwt[3]);
				}
			} /*end interpolate this point*/
		} /*end only do this point if there's a DEM value*/
	} /*end for all points to be interpolated*/
	
	printf("  Done!\n\n");
	
	/*WRITE DEM and RDR elevations to output, if both have values!*/
	/*X,Y,Z ASCII Lists**********************************************/
	DEMfile = fopen(Filenames[2], "w");
	RDRfile = fopen(Filenames[3], "w");
	
	if (DEMfile == NULL) {
		printf("Cannot open DEM output file=[%s]:[%s]! Exiting.\n",
		Filenames[2], strerror(errno));
		return(-1);
	}
	if (RDRfile == NULL) {
		printf("Cannot open Radar output file=[%s]:[%s]! Exiting.\n",
		Filenames[3], strerror(errno));
		return(-1);
	}

	j=0;
	for (i=0;i<DEMinboundsct;i++){
		if (RDRinterp[i].elevation!=-9999) {
			j++;
			fprintf(DEMfile, "%0.3f\t%0.3f\t%0.6f\n", 
				               cropDEMvals[i].easting,
				               cropDEMvals[i].northing,
				               cropDEMvals[i].elevation);
			fprintf(RDRfile, "%0.3f\t%0.3f\t%0.6f\n", 
				               RDRinterp[i].easting,
				               RDRinterp[i].northing,
				               RDRinterp[i].elevation);
				               
		}
	}

	fclose(RDRfile);
	fclose(DEMfile);
	printf("%u points written to ASCII Output files:\n  %s (DEM)\n  %s (RDR)\n",
	       j, Filenames[2],Filenames[3]);
	
	/*Calculate simulation time elapsed, and print it.*/
	endTime = time(NULL);
	printf("\nElapsed Time approximately %u seconds.\n\n",
	       (unsigned)(endTime - startTime));
	
	return(0);
}
