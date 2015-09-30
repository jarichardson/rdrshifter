#include "structs.h"

double *DEM_LOADER(char *DEMfilename, point **dempts) {
/*Module DEM_LOADER_GDAL
	Accepts a file name and a null data grid
	Checks for validity of raster file
	loads raster file metadata into array (DEMGeoTransform)
	Loads Raster Data into Global Data Grid depending on Raster Type given:
		TOPOG: Digital Elevation Model, put into DEMgrid.dem_elev
		T_UNC: DEM Uncertainty        , put into DEMgrid.elev_uncert
		RESID: Flow Residual Thickness, put into DEMgrid.residual
	
	passes data grid back through a pointer
	returns metadata array, or a NULL Double for errors
*/
	
	GDALDatasetH    DEMDataset;/*The raster file*/
	GDALDriverH     DEMDriver; /*essentially the raster File Type*/
	GDALRasterBandH DEMBand;
	double          *DEMGeoTransform;
	float           *pafScanline;
	point           *DEMgrid;

	unsigned        YOff;
	unsigned        i,j;
	double          easting,northing;
	
	
	if((DEMGeoTransform = malloc (sizeof (double) * 6))==NULL) {
		printf("ERROR [DEM_LOADER]: Out of Memory creating Metadata Array!\n");
		return((double*) -1);
	}
	
	/*DEMGeoTransform[0] lower left x
	  DEMGeoTransform[1] w-e pixel resolution
	  DEMGeoTransform[2] number of cols, assigned manually in this module 
	  DEMGeoTransform[3] lower left y
	  DEMGeoTransform[4] number of lines, assigned manually in this module
	  DEMGeoTransform[5] n-s pixel resolution (negative value) */
	
	GDALAllRegister();
	
	/*Open DEM raster file*/
	DEMDataset = GDALOpen( DEMfilename, GA_ReadOnly );
	if(DEMDataset==NULL){
		printf("ERROR [DEM_LOADER]: DEM file could not be loaded!\n");
		printf("  File:  %s\n", DEMfilename);
		return((double*) NULL);
	}
	
	/*Make sure Projection metadata is valid (that the raster is a valid raster)*/
	if( GDALGetProjectionRef( DEMDataset ) == NULL )
	{
		printf("ERROR [DEM_LOADER]: DEM Projection metadata could not be read!\n");
		printf("  File:  %s\n", DEMfilename);
		return((double*) NULL);
	}
	
	/*Read DEM raster metadata into DEMGeoTransform*/
	if( GDALGetGeoTransform( DEMDataset, DEMGeoTransform ) != CE_None )
	{
		printf("ERROR [DEM_LOADER]: DEM GeoTransform metadata could not be read!\n");
		printf("  File:  %s\n", DEMfilename);
		return((double*) NULL);
	}
	
	/*Find out the raster type (e.g. tiff, netcdf)*/
	DEMDriver = GDALGetDatasetDriver( DEMDataset );
	
	/*Modify Metadata for use in this program*/
	DEMGeoTransform[5] = -1 * DEMGeoTransform[5]; /*row height*/
	/*DEMGeoTransform[1]; col width*/
	DEMGeoTransform[4] = GDALGetRasterYSize( DEMDataset ); /*number of rows*/
	DEMGeoTransform[2] = GDALGetRasterXSize( DEMDataset ); /*number of cols*/
	DEMGeoTransform[3] -= (DEMGeoTransform[5] * DEMGeoTransform[4]);/*Lower left Y value*/
	/*DEMGeoTransform[0]; Lower left X value*/
	
	/*Output DEM metadata to the screen*/
	printf("\nRaster read from %s file successfully.\n\nRaster Information:\n",
	       GDALGetDriverLongName( DEMDriver ) );
	printf("  File:              %s\n", DEMfilename);
	printf("  Lower Left Origin: (%.6f,%.6f)\n",
	       DEMGeoTransform[0], DEMGeoTransform[3]);
	printf("  GMT Range Code:    -R%.3f/%.3f/%.3f/%.3f\n",
	        DEMGeoTransform[0],
	       (DEMGeoTransform[0]+((DEMGeoTransform[2]-1)*DEMGeoTransform[1])),
	        DEMGeoTransform[3],
	       (DEMGeoTransform[3]+((DEMGeoTransform[4]-1)*DEMGeoTransform[5])));
	printf("  Pixel Size:        (%.6f,%.6f)\n",
	       DEMGeoTransform[5], DEMGeoTransform[1]);
	printf("  Grid Size:         (%u,%u)\n",
	       (unsigned) DEMGeoTransform[4], (unsigned) DEMGeoTransform[2]);
	
	/*Create 2D global data grid from DEM information*/
	printf("\nAllocating Memory for Point list (%u values)...\n",
	       ((unsigned) DEMGeoTransform[4] * (unsigned) DEMGeoTransform[2]));
	
	/*Load elevation information into DEMBand. DEM should be only band in raster*/
	DEMBand = GDALGetRasterBand(DEMDataset, 1);
	/*Allocate memory the length of a single row in the DEM*/
	pafScanline = (float *) CPLMalloc(sizeof(float)*DEMGeoTransform[2]);
	
	/*allocate memory for data grid*/
	if((*dempts = malloc (sizeof (point) * ((unsigned) DEMGeoTransform[4] * (unsigned) DEMGeoTransform[2])))==NULL) {
		printf("ERROR [DEM_LOADER]: Out of Memory creating Metadata Array!\n");
		return((double*) -1);
	}
		
	DEMgrid = *dempts; /*Assign pointer position to DEMgrid variable*/
	
	for(i=0;i<DEMGeoTransform[4];i++) { /*For each row*/
		/*calculate row so bottom row is read into data array first*/
		YOff = (DEMGeoTransform[4]-1) - i;

		/*Read elevation data from row in input raster*/
		if((GDALRasterIO(DEMBand, GF_Read, 0, YOff, DEMGeoTransform[2], 1, 
				             pafScanline, DEMGeoTransform[2], 1, GDT_Float32, 0, 0))
			 != CE_None) {
			printf("\nERROR [DEM_LOADER]: DEM Elevation Data could not be read!\n");
			printf("  File:  %s\n", DEMfilename);
			return((double*) NULL);
		}
		
		
		/*Write elevation data column by column into list array using DEMgrid variable and X,Y values*/
		/*Calculate Y value, i is the current row (max=#cols-1, heighest row; min=0, lowest row).*/
		northing = DEMGeoTransform[3] + (DEMGeoTransform[5] * i);
		
		for(j=0;j<DEMGeoTransform[2];j++) {
			/*Calculate X value*/
			easting = DEMGeoTransform[0] + (DEMGeoTransform[1] * j);
			DEMgrid[(unsigned) ((i*DEMGeoTransform[2])+j)].easting   = easting;
			DEMgrid[(unsigned) ((i*DEMGeoTransform[2])+j)].northing  = northing;
			DEMgrid[(unsigned) ((i*DEMGeoTransform[2])+j)].elevation = pafScanline[j];
		}
	}
	
	printf("\n                  Elevation Model Loaded into list.\n\n");
	
	CPLFree(pafScanline);
	return(DEMGeoTransform);

}
