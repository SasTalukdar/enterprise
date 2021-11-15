Star(t)Date 2021.07.01.12.30  
# Entry 1 : To install enterprise use "pip3 install ."  
-----------------------------------------------------------------------------------------------  
-----------------------------------------------------------------------------------------------  
## Captains log: DeckOne is used for various plots  
         
###        SpatialPlots : can be used to initialize spatial plot axes and domain
               
####            initialize(lat1,lat2,lon1,lon2,res=0,lb_size=14,rotate=0,linewidth=0.8,dec_lim=0,shape_file_path=None,title=None):
                        
                Return: geoAxis
                        
                        # Initialize a geoaxis with background map for ploting
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

####            initialize_subplots(lat1, lat2, lon1, lon2, plot_no=1, res=0, lb_size=14, rotate=0, linewidth=0.8, dec_lim=0, shape_file_path=None, add_label=False, titles=None,
fig_sg=None,x_nos=None,y_nos=None, coastlines=False):

Return: figure, list of geoAxis

      ** Initialize subplots **
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

####           polar_plot(lon,lat,val,lon0=None,lat0=None,r_lim=None,t_lim=None,th_off=np.pi/2.0,
####                      t_tick=[0,np.pi*0.5,np.pi,np.pi*1.5],r_tick=None,r_pos=None,c_map='jet')

               Return: figure, Axis

                        # Polar contour plot using lat, lon matrices
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
