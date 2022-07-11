<img src="https://github.com/Doktor-Hu/RACMO_Tool/blob/main/img/RACMOTOOL.png" alt="drawing" width="200"/> 

# 1. Introduction

RACMO Tool is a Python-and GDAL-based package for pre-processing, read, plot standard RACMO2 NetCDF files produced by IMAU UU. It is sepecifically designed for RACMO2 users from non-modeller domain, e.g., Remotes Sensing, geo-statistics. Or to make it easier to complie to QGIS.

#### Warning: GDAL version >=3.1.4 is required, checking your GDAL version!

# 2. The functionalities includes:

1) Convert NetCDF file in to GeoTiff with a default EPSG:3031 projection.

   -- fn: Input filename
  
   -- param: Parameter short name (see file name or pdf);
   
   -- res: Input file resolution in [m] 27000 or 5500;
   
   -- out_fn_var: Output filename;
   
   -- exe: Boolean, if the code short be executed or only print the commandline.
```
import RACMOTool as R2Tool
R2Tool.NC_to_TIFF(fn, param, res, out_fn_var, exe=False)
```
2) Extracting time series from a point(s) ot shapefile.
3) Batch processing

```
RACMOTool.py ..../Dir 27000 clcov.KNMI-2001.ANT27.ERAINx_RACMO2.3p2.DD.nc N
```

The code are modified from ... and ..., and automated by Dr. Zhongyang Hu from IMAU/UU.
