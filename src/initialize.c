#include "structs.h"

int INITIALIZE(char *CFGfilename, char ***Filenames, double *or_east, double *or_north, 
               double *or_elev, double *tr_east, double *tr_north, double *tr_elev, 
               double *ro_east, double *ro_north, double *ro_elev, double *param_scale, 
               double *param_DEMnodata, double *param_RDRnodata) {
/*Module INITIALIZE
	Accepts a configuration file and returns model variables:
	
	string DEM_RASTER
	string RDR_RASTER
	string OUTPUT_DEM 
	string OUTPUT_RDR 

	double ORIGIN_X
	double ORIGIN_Y
	double ORIGIN_Z

	double TRANSLATE_X
	double TRANSLATE_Y
	double TRANSLATE_Z
	double ROTATE_X
	double ROTATE_Y
	double ROTATE_Z
	double SCALE
	
	Checks at end for configuration file errors, where mandatory parameters were
	  not assigned.
*/

	unsigned maxLineLength = 256;
	char     line[256];             /*Line string from file       */
	char     var[64];               /*Parameter Name  in each line*/
	char     value[256];            /*Parameter Value in each line*/
	int      i;
	char     *ptr;
	char     **INITFiles;      /*working array for filenames in this module     */
	unsigned INITFilesCount = 4;   /*number of possible input files             */
	FILE     *Opener;     /*Dummy File variable to test valid output file paths */
	unsigned outputs = 0; /*Number of output types called for (must have 2)   */
	unsigned DOFs = 0;    /*Number of DOF to work with (must have >=1)*/
	unsigned origins = 0;
	FILE *ConfigFile;
	
	
	/*Allocate filename pointers*/
	if((*Filenames=(char**)malloc((unsigned)(INITFilesCount+1)*sizeof(char*)))==NULL){
		printf("ERROR [INITIALIZE]:\n");
		printf("   NO MORE MEMORY: Tried to allocate memory for 4 filenames\n");
		return(-1);
	}
	INITFiles = *Filenames;
	for(i=0;i<(INITFilesCount+1);i++){
		if((INITFiles[i]=(char*)malloc(sizeof(char)*(maxLineLength+1)))==NULL) {
			printf("\n[INITIALIZE] Out of Memory assigning filenames!\n");
			return(-1);
		}
	}
	
	printf("Reading in Parameters...\n");
	
	/*open configuration file*/
	ConfigFile = fopen(CFGfilename, "r");
	if (ConfigFile == NULL) {
		printf("\nERROR [INITIALIZE]: Cannot open configuration file=[%s]:[%s]!\n",
		       CFGfilename,
		       strerror(errno));
		return(-1);
	}

	/* use each line to compare to needed values*/
	while (fgets(line, maxLineLength, ConfigFile) != NULL) {
		/*if first character is comment, new line, space, return to next line*/
		if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') continue;
		
		/*print incoming parameter*/
		sscanf (line,"%s = %s",var,value); /*split line into before ' = ' and after*/
		printf("%20s = %-33s ",var,value); /*print incoming parameter value*/
		
		/*INPUT FILES AND GLOBAL MODEL PARAMETERS**********************************/
		/*INPUT DEM FILE*/
		if (!strncmp(var, "DEM_RASTER", strlen("DEM_RASTER"))) {
			strcpy(INITFiles[0],value);
			printf("[assigned]\n");
		}
		
		/*INPUT RADAR DEM FILE*/
		else if (!strncmp(var, "RDR_RASTER", strlen("RDR_RASTER"))) {
			strcpy(INITFiles[1],value);
			printf("[assigned]\n");
		}
		
		/*OUTPUT FILES*************************************************************/
		/*OUTPUT ASCII X,Y,Z DEM FILE*/
		else if (!strncmp(var, "OUTPUT_DEM", strlen("OUTPUT_DEM"))) {
			strcpy(INITFiles[2],value);
			
			Opener = fopen(INITFiles[2], "w");
			if (Opener == NULL) {
				printf("\nERROR [INITIALIZE]: Failed to create an output file at [%s]:[%s]!\n",
				       INITFiles[2],
				       strerror(errno));
				return(-1);
			}
			else {
				(void) fclose(Opener);
				outputs++;
				printf("[assigned]\n");
			}
		}
		
		/*MATCHING OUTPUT ASCII X,Y,Z RADAR FILE*/
		else if (!strncmp(var, "OUTPUT_RDR", strlen("OUTPUT_RDR"))) {
			strcpy(INITFiles[3],value);
			
			Opener = fopen(INITFiles[3], "w");
			if (Opener == NULL) {
				printf("\nERROR [INITIALIZE]: Failed to create an output file at [%s]:[%s]!\n",
				       INITFiles[3],
				       strerror(errno));
				return(-1);
			}
			else {
				(void) fclose(Opener);
				outputs++;
				printf("[assigned]\n");
			}
		}
		
		/*NoDATA*/
		else if (!strncmp(var, "RDR_NODATA", strlen("RDR_NODATA"))) {
				*param_RDRnodata = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Easting Translation (RDR_NODATA) is not a number!\n");
					return(-1);
				}
				printf("[assigned]\n");
		}
		
		else if (!strncmp(var, "DEM_NODATA", strlen("DEM_NODATA"))) {
				*param_DEMnodata = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Easting Translation (DEM_NODATA) is not a number!\n");
					return(-1);
				}
				printf("[assigned]\n");
		}
		
		/*6-DOF SHIFT PARAMETERS**************************************************/
		/*ORIGIN LOCATION*/
		else if (!strncmp(var, "ORIGIN_X", strlen("ORIGIN_X"))) {
				*or_east = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Origin Easting (ORIGIN_X) is not a number!\n");
					return(-1);
				}
				origins+=1;
				printf("[assigned]\n");
		}
		else if (!strncmp(var, "ORIGIN_Y", strlen("ORIGIN_Y"))) {
				*or_north = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Origin Northing (ORIGIN_Y) is not a number!\n");
					return(-1);
				}
				origins+=1;
				printf("[assigned]\n");
		}
		else if (!strncmp(var, "ORIGIN_Z", strlen("ORIGIN_Z"))) {
				*or_elev = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Origin Elevation (ORIGIN_Z) is not a number!\n");
					return(-1);
				}
				origins+=1;
				printf("[assigned]\n");
		}
		
		/*Translation*/
		else if (!strncmp(var, "TRANSLATE_X", strlen("TRANSLATE_X"))) {
				*tr_east = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Easting Translation (TRANSLATE_X) is not a number!\n");
					return(-1);
				}
				DOFs+=1;
				printf("[assigned]\n");
		}
		else if (!strncmp(var, "TRANSLATE_Y", strlen("TRANSLATE_Y"))) {
				*tr_north = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Northing Translation (TRANSLATE_Y) is not a number!\n");
					return(-1);
				}
				DOFs+=1;
				printf("[assigned]\n");
		}
		else if (!strncmp(var, "TRANSLATE_Z", strlen("TRANSLATE_Z"))) {
				*tr_elev = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Elevation Translation (TRANSLATE_Z) is not a number!\n");
					return(-1);
				}
				DOFs+=1;
				printf("[assigned]\n");
		}
		
		/*Rotation*/
		else if (!strncmp(var, "ROTATE_X", strlen("ROTATE_X"))) {
				*ro_east = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Rotation about the X-axis (ROTATE_X) is not a number!\n");
					return(-1);
				}
				DOFs+=1;
				printf("[assigned]\n");
		}
		else if (!strncmp(var, "ROTATE_Y", strlen("ROTATE_Y"))) {
				*ro_north = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Rotation about the Y-axis (ROTATE_Y) is not a number!\n");
					return(-1);
				}
				DOFs+=1;
				printf("[assigned]\n");
		}
		else if (!strncmp(var, "ROTATE_Z", strlen("ROTATE_Z"))) {
				*ro_elev = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Rotation about the Z-axis (ROTATE_Z) is not a number!\n");
					return(-1);
				}
				DOFs+=1;
				printf("[assigned]\n");
		}
		
		/*Scale*/
		else if (!strncmp(var, "SCALE", strlen("SCALE"))) {
				*param_scale = strtod(value,&ptr);
				if (ptr == value) { /*NOT A NUMBER*/
					printf("\nERROR [INITIALIZE]: Scale Factor (SCALE) is not a number!\n");
					return(-1);
				}
				DOFs+=1;
				printf("[assigned]\n");
		}
		
		else {
			printf("[not assigned]\n");
			continue;
		}
	}
	
	
	/*Check for missing parameters*/
	if(!*param_scale) { /*Elevation uncertainty is either missing or is 0.*/
		*param_scale = 1.0;
		printf("SCALE = 1.0: Radar values will not be spatially scaled.\n");
	}
	if(!strcmp(INITFiles[0],"")) { /*DEM Filename is missing.*/
		printf("\nERROR [INITIALIZE]: No DEM Filename Given!!\n");
		return(-1);
	}
	if(!strcmp(INITFiles[1],"")) { /*Radar DEM Filename is missing.*/
		printf("\nERROR [INITIALIZE]: No Radar DEM Filename Given!!\n");
		return(-1);
	}
	if(outputs!=2) { /*No Output Filenames are given.*/
		printf("\nERROR [INITIALIZE]: Not enough Output Filenames Given!!\n");
		return(-1);
	}
	if(origins!=3) { /*No Output Filenames are given.*/
		printf("\nERROR [INITIALIZE]: No Origin location given!!\n");
		return(-1);
	}
	if(!DOFs) { /*No Output Filenames are given.*/
		printf("\nERROR [INITIALIZE]: No Rotation/Translation/Scale Modifiers given!!\n  At least use 1.\n");
		return(-1);
	}
	
	
	(void) fclose(ConfigFile);
	return(0);
}
