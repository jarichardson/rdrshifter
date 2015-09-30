#ifndef _DOF_STRUCTS_
#define _DOF_STRUCTS_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <float.h>
#include <string.h>
#include "structs.h"  /* Global Structures and Variables*/
#include "gdal.h"     /* GDAL */
#include "cpl_conv.h" /* GDAL for CPLMalloc() */

/*3D Points*/
typedef struct point {
	double northing;
	double easting;
	double elevation;
} point;

/*Global Variables*/
time_t startTime;
time_t endTime;

int INITIALIZE(char*,char***,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

double *DEM_LOADER(char*, point**);

/*      ,``                 .`,   
      `,.+,               ``+`    
   ,.`#++``               ,;++`,  
 , +++'`.                  ,++++ `
` +++++++;.``  ,   ,     `#+++++ `
``+++++++++++`,   ,'++++++++++'`` 
  .,`  `++++ .   .++++++;.  ,,    
     ,+++++.   ` +++++ ,          
     ,`+++;, ``'++++`,            
    ,;+++`, . +++++,              
  . ++++;;:::++++',               
   ,'++++++++++#`,                
    ,.....````,`                  
                                  
        GO BULLS                */

#endif /*#ifndef _DOF_STRUCTS_*/
