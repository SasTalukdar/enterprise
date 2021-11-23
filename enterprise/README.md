Star(t)Date 2021.07.01.12.30 
## Captains log: DeckOne is used for various plots  
         
###        SpatialPlots : can be used to initialize spatial plot axes and domain
               
####   initialize(lat1, lat2, lon1, lon2, res=0, lb_size=14, rotate=0, linewidth=0.8, dec_lim=0, shape_file_path=None, title=None):

Return: geoAxis

Initialize a geoaxis with background map for ploting
      
      lat1 (float): minimum latitude
      lat2 (float): maximum latitude
      lon1 (float): minimum longitude
      lon2 (float): maximum longitude
      res  (float): (optional) distance between two consecutive of x_ticks/y_ticks
      lb_size (float): (optional) label size
      rotate (float): (optional) rotation angle for x and y tick labels
      linewidth (float): (optional) linewidth
      dec_lim (integer): (optional) decimal places to which ticks are rounded off
      shape_file_path (string/geometry):  shape file to be plotted in the figure (can be path or a geometry object)
      title (string): (optional) title for the figure

####   initialize_subplots(lat1, lat2, lon1, lon2, plot_no=1, res=0, lb_size=14, rotate=0, linewidth=0.8, dec_lim=0, shape_file_path=None, add_label=False, titles=None, fig_sg=None,x_nos=None,y_nos=None, coastlines=False):

Return: figure, list of geoAxis

Initialize subplots 

      lat1 (float): minimum latitude
      lat2 (float): (float) maximum latitude
      lon1 (float): minimum longitude
      lon2 (float): maximum longitude
      plot_no (integer): total number of subplots
      res (float): (optional) distance between two consecutive of x_ticks/y_ticks
      lb_size (float): (optional) label size
      rotate (float): (optional) rotation angle for x and y tick labels
      linewidth (float): (optional) linewidth
      dec_lim (integer): (optional) decimal places to which ticks are rounded off
      shape_file_path (string/geometry): shape file to be plotted in the figure (can be path or a geometry object)
      add_label (True/False): (optional) label for each subplot
      titles (list of string): (optional) titles for each plot
      fig_sg (tuple of float): (optional) figure size for the plot
      x_nos (integer): (optional) number of plots in the x direction
      y_nos (integer): (optional) number of plots in the y direction
      coastlines (True/False): (optional) add coastlines

####   polar_plot(lon, lat, val,lon0=None, lat0=None, r_lim=None, t_lim=None, th_off=np.pi/2.0, t_tick=[0,np.pi*0.5,np.pi,np.pi*1.5], r_tick=None, r_pos=None, c_map='jet')

Return: figure, Axis

Polar contour plot using lat, lon matrices

      lon (array): longitude array
      lat (array): latitude array
      val (2d array): values to be ploted
      lon0 (float): (optional) center longitude
      lat0 (float): (optional) center latitude
      r_lim (float): (optiona) maximum limit for radius
      t_lim (float): (optional) maximum limit for theta in radians
      th_off (float): (optional) theta offset in radian
      t_tick (array): (optional) theta ticks
      r_tick (array): (optional) radial ticks
      r_pos (float): (optional) position of r ticks in radian
      c_map (string): (optional) colormap for plotting

###        AnalyticalPlots : can be used for analytical plots
####  ContouredScatter(x, y, bins=100, level=30, cmap='jet', savename='contoured_scatter', limits=None)

Return None

Contoured scatter plot for large sample sizes

         x (array): array
         y (array): array
         bins (integer): bins
         level (integer): number of contour levels
         camp (string): colormap
         savename (string): filename to be saved as
         limits (array of size 4): the limits for the plots ([x_min, x_max, y_min, y_max])


## Captains log: DeckTwo contains WRF related functions 

###        wps : This class contains scripts related to WPS

####   namelist_plot(namelist_path)

Return: geoAxis

Plot the domain set up using the namelist.wps file

      namelist_path (string): path of the namelist file

#### inject_lulc(input_lulc,exchange_pairs_arr,output_lulc='geo_em.d01.nc')

Return: 3D array as a replacement for the 'LU_INDEX' variable in the wrf input file

Replace the default LULC values with updated values

	input_lulc (GeoTiff (.tiff) file): updated LULC values for the domain
	exchange_pairs_arr (list of tuples): an array of to be interchanged values as an array of tuples.
	output_lulc (string): (optional) name of the output file to be created

###	  assimilation : This class contains scripts related to data assimilation

####	insat2littler(data_paths)

Return: file

Create observation file in the little-R format using provided INSAT profile data

	data_paths (list of string): list of INSAT data files

####	show_assimilated_points(diag_path,namelist_path=None,res=0,s=np.pi/30,alpha=0.7,lb_size=11)

Return: geoAxis

Show the assimilated point locations in the domain

	diag_path (string): path to the diagnostic file
	namelist_path (string): (optional) path to the namelist.input file for the consideration of domain boundary
	res (float): (optional) the distance between two respective x/y tick
	s (float): (optiona) the size of the location markers
	alpha (float): (optional) the transparency of the location markers
	lb_size (float): label size

## Captains log: DeckThree is used for analytical fuctions

###	interpolate : This class contains interpolation related functions

####	regrid(data, out_x, out_y)

Return: 2-dimentional array

Interpolate a given matrix into any shape or size

	data (2D array): input data
	out_x (integer): number of grids in the x axis of the interpolated array
	out_y (integer): number of grids in the y axis of the interpolated array

### statistics : This class contains statistical fuctions

####	corr(obs,exp)

Return: [float, float]

Calculate Pearsonr correlation coefficient between two 2-D data arrays. Returns [correlation, p-value]

	obs (2D array): observation data
	exp (2D array): experimental data

####    rmse(obs,exp)

Return: float

Calculate root mean square error between two 2-D data arrays

        obs (2D array): observation data
        exp (2D array): experimental data

####    bias(obs,exp)

Return: float

Calculate bias between two 2-D data arrays

        obs (2D array): observation data
        exp (2D array): experimental data

####	con_int(data, confidence = 0.95)

Return: float

Calculate confidence interval for a given data array

	data (1D array): data
	confidence (float): confidence level

####	remove_outlier(data)

Return: array

Replace outlier values from a given data array with nan using IQR method
following http://mathcenter.oxford.emory.edu/site/math117/shapeCenterAndSpread/

	data (1D array): data

### cyclone : This class contains cyclone related analysis

####	vel2polvel(u,v,lat,lon,lat0=None,lon0=None)

Return: (array, array)

Convert zonal and meridional wind to radial and azimuthal wind (v_r,v_th)

	u (2D array): meridional wind data array
	v (2D array): zonal wind data array
	lat (array): latitude data array
	lon (array): longitude data array
	lat0 (float): (optional) central latitude
	lon0 (float): (optional) central longitude
	If lat0, lon0 are not provided, the center is calculated based on minimum wind speed criteria

####	agm_avg(datapath,t=0)

Return: None

Plot azimuthal average of radial and azimuthal wind speed from a given wrf output file (.nc)
	datapath (string): path to the wrf output file
	t (integer): (optional) time step to be plotted


# The contruction of the documentation is in progress. Please visit again in a few days.
