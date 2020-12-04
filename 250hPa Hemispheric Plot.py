from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
from netCDF4 import num2date
import numpy as np
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS

def queryData(ncssObj, dataType):
	# Create NCSS Query for heights
	ncssQuery = ncss.query().time(now).accept('netcdf4')
	ncssQuery.lonlat_box(0, 360, 0, 90)# Set the lat/lon box for the data you want to pull in.
	ncssQuery.vertical_level(25000)# Set desired level 50000 = 50000 Pa = 500 hPa

	if dataType == "height":
		ncssQuery.variables('Geopotential_height_isobaric').add_lonlat()
	elif dataType == "wind":
		ncssQuery.variables('u-component_of_wind_isobaric',
                       'v-component_of_wind_isobaric').add_lonlat()

	#gfsdata_hght query is var=Geopotential_height_isobaric&time=2020-11-28T21%3A47%3A58.499031&west=0&east=360&south=0&north=90&accept=netcdf4&addLatLon=True&vertCoord=25000

	# Get the data using the query
	data = ncss.get_data(ncssQuery)
	return data

def plotAll(Z_250, wspd250):
	datacrs = ccrs.PlateCarree()
	plotcrs = ccrs.NorthPolarStereo(central_longitude=-100.0)

	# Make a grid of lat/lon values to use for plotting with Basemap.
	lons, lats = np.meshgrid(longitudes, latitudes)

	fig = plt.figure(1, figsize=(12., 13.))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02],
                       bottom=.07, top=.99, hspace=0.01, wspace=0.01)

	axes = plt.subplot(gs[0], projection=plotcrs)
	axes.set_title('250-hPa Geopotential Heights (m)', loc='left')
	#ax.set_title('VALID: {}'.format(vtimes[0]), loc='right')

	#   ax.set_extent([west long, east long, south lat, north lat])
	axes.set_extent([-180, 180, 10, 90], datacrs)
	axes.coastlines('50m', edgecolor='black', linewidth=0.5)
	axes.add_feature(cfeature.STATES, linewidth=0.5)

	contourLvlHeights = np.arange(9000, 12000, 120) #start at 9000m, end at 12000m, 120m interval
	contours = axes.contour(lons, lats, Z_250, contourLvlHeights, colors='k',
    	            linewidths=1.0, linestyles='solid', transform=datacrs)
	plt.clabel(contours, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

	contourLvlWinds = np.arange(50, 230, 20) #start at 50 knots, end at 230, 20 knot interval
	colormap = plt.cm.get_cmap('BuPu')
	contourfill = axes.contourf(lons, lats, wspd250, contourLvlWinds, cmap=colormap, transform=datacrs)
	cax = plt.subplot(gs[1])
	cbar = plt.colorbar(contourfill, cax=cax, orientation='horizontal', extend='max', extendrect=True)

	plt.show()

def setupPlot():
	
	fig = plt.figure(1, figsize=(12., 13.))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02],
                       bottom=.07, top=.99, hspace=0.01, wspace=0.01)
	axes = plt.subplot(gs[0], projection=ccrs.NorthPolarStereo(central_longitude=-100.0))
	axes.set_title('250-hPa Geopotential Heights (m)', loc='left')

	#   ax.set_extent([west long, east long, south lat, north lat])
	axes.set_extent([-180, 180, 10, 90], ccrs.PlateCarree())
	axes.coastlines('50m', edgecolor='black', linewidth=0.5)
	axes.add_feature(cfeature.STATES, linewidth=0.5)

	return axes, gs

def setupContourLvlHeights(lons, lats, axes, Z_250, plt):
	contourLvlHeights = np.arange(9000, 12000, 120) #start at 9000m, end at 12000m, 120m interval
	contours = axes.contour(lons, lats, Z_250, contourLvlHeights, colors='k',
    	            linewidths=1.0, linestyles='solid', transform=ccrs.PlateCarree())
	plt.clabel(contours, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

def setupContourLvlWinds(plt, axes, lons, lats, wspd250):
	contourLvlWinds = np.arange(50, 230, 20) #start at 50 knots, end at 230, 20 knot interval
	colormap = plt.cm.get_cmap('BuPu')
	contourfill = axes.contourf(lons, lats, wspd250, contourLvlWinds, cmap=colormap, transform=ccrs.PlateCarree())
	cax = plt.subplot(gs[1])
	cbar = plt.colorbar(contourfill, cax=cax, orientation='horizontal', extend='max', extendrect=True)


# Latest GFS Dataset
cat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/'
                 'NCEP/GFS/Global_0p5deg/latest.xml')
#cat.datasets = [GFS_Global_0p5deg_20201128_1200.grib2]
best_ds = list(cat.datasets.values())[0]
ncss = NCSS(best_ds.access_urls['NetcdfSubset']) #NetCDF subset service object

now = datetime.utcnow()

data_hght = queryData(ncss, "height")
data_wind = queryData(ncss, "wind")

latitudes = data_hght.variables['lat'][:]
longitudes = data_hght.variables['lon'][:]
heights250hPa = data_hght.variables['Geopotential_height_isobaric'][:]

# Smooth the 250-hPa heights using a gaussian filter from scipy.ndimage
hgt_250, longitudes = cutil.add_cyclic_point(heights250hPa,
                                      coord=longitudes)
Z_250 = ndimage.gaussian_filter(hgt_250[0, 0, :, :], sigma=3, order=0)

u250 = (units(data_wind.variables['u-component_of_wind_isobaric'].units) *
        data_wind.variables['u-component_of_wind_isobaric'][0, 0, :, :])
v250 = (units(data_wind.variables['v-component_of_wind_isobaric'].units) *
        data_wind.variables['v-component_of_wind_isobaric'][0, 0, :, :])

u250 = u250.units * cutil.add_cyclic_point(u250)
v250 = v250.units * cutil.add_cyclic_point(v250)
wspd250 = mpcalc.wind_speed(u250, v250).to('knots')

# Make a grid of lat/lon values to use for plotting with Basemap.
lons, lats = np.meshgrid(longitudes, latitudes)

#plotAll(Z_250, wspd250)
axes, gs = setupPlot()
setupContourLvlHeights(lons, lats, axes, Z_250, plt)
setupContourLvlWinds(plt, axes, lons, lats, wspd250)
plt.show()

