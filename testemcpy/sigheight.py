import numpy as np
import matplotlib.pyplot as plt
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded
import netCDF4 as nc
import argparse
import os

def gen_figure(date, outpath):
    variable='HTSGW_surface'
    #Replace fill values with nan, returns corrected values, latitude, longitude
    def correct(dataset, variable, lats, lons, start=0, stop=1):
        fill_val = dataset.variables[variable].__dict__['_FillValue']
        time = dataset.variables['time'][start:stop]
        original = np.array(dataset.variables[variable][start:stop])
        corrected = np.full(original.shape, np.nan)
        for i in range(len(time)):
            for j in range(len(lats)):
                for k in range(len(lons)):
                    if original[i][j][k] == fill_val:
                        continue
                    else:
                        corrected[i][j][k] = original[i][j][k]
        return corrected

        
    def get_wind(model):
        model_sigheight=np.full((129, 721, 1440), np.nan)
        
        for i in range(129):
            hour=3*i
            if len(str(hour))==1:
                i3="00"+str(hour)
            elif len(str(hour))==2:
                i3="0"+str(hour)
            else:
                i3=str(hour)
            f = "/work2/noaa/marine/ifox/{dtype}/{dtype}{date}00/COMROOT/{dtype}{date}00/gfs.{date}/00/products/wave/gridded/gfswave.t00z.global.0p25.f{hour}.grib2.nc".format(dtype=model, hour=i3, date=date)
            dataset=nc.Dataset(f)
            modlats = np.array(dataset.variables['latitude'][:])
            modlons = np.array(dataset.variables['longitude'][:])

            sigheight=correct(dataset, variable, modlats, modlons)
            
            model_sigheight[i]=sigheight[0]
        return model_sigheight, modlats, modlons
            
    coupled_height, lats, lons=get_wind('onewaycpl')
    control_height, lats, lons=get_wind('control')
    
    coupled_control_diff=coupled_height-control_height
    '''
    print(np.nanmin(analysis_wind_speed), np.nanpercentile(analysis_wind_speed, 1), np.nanpercentile(analysis_wind_speed, 99), np.nanmax(analysis_wind_speed))
    print(np.nanmin(model_wind_speed), np.nanpercentile(model_wind_speed, 1), np.nanpercentile(model_wind_speed, 99), np.nanmax(model_wind_speed))

    print(np.nanmin(analysis_u), np.nanpercentile(analysis_u, 1),np.nanpercentile(analysis_u, 99), np.nanmax(analysis_u))
    print(np.nanmin(model_u), np.nanpercentile(model_u, 1), np.nanpercentile(model_u, 99), np.nanmax(model_u))
    '''
        
    X, Y=np.meshgrid(lats, lons)

    def map(data, minimum, maximum, cmap, name):

        # Create gridded map object
        gridded = MapGridded(X, Y, np.transpose(data))
        gridded.cmap = cmap
        gridded.vmin = minimum
        gridded.vmax = maximum

        # Create plot object and add features
        plot = CreatePlot()
        plot.plot_layers = [gridded]
        plot.projection = 'plcarr'
        plot.domain = 'global'
        plot.add_map_features(['coastline'])
        #plot.add_xlabel(xlabel='longitude')
        #plot.add_ylabel(ylabel='latitude')
        plot.add_title(label=name, loc='center', fontsize=14)
        plot.add_grid()
        plot.add_colorbar(label='m', fontsize=12, extend='neither')
        return plot

    def fig(data, name):
            fig1 = CreateFigure()
            fig1.plot_list = [data]
            fig1.create_figure()
            fig1.save_figure(name)
            fig1.close_figure()


    datamin=np.nanmin((coupled_height, control_height))
    datamax=np.nanmax((coupled_height, control_height))
    diffmax=np.nanmax(np.absolute(coupled_control_diff))
    print(diffmax, np.nanpercentile(np.absolute(coupled_control_diff), 99))

    for i in range(129):
        fhour=3*i
        if len(str(fhour))==1:
            fhour3='00'+str(fhour)
        elif len(str(fhour))==2:
            fhour3='0'+str(fhour)
        else:
            fhour3=fhour
        for j in [6]:
            coupled_plot=map(coupled_height[i], 0, 15 , 'turbo', 'One-way coupled forecast hour {fhour} significant wave height'.format(fhour=fhour3))
            control_plot=map(control_height[i], 0, 15, 'turbo', 'Two-way coupled forecast hour {fhour} significant wave height'.format(fhour=fhour3))
            coupled_control_plot=map(coupled_control_diff[i], -1*j, j, 'seismic', 'One-way coupled - two-way coupled forecast hour {fhour} significant wave height'.format(date=date, fhour=fhour3))
                                    
            fig(control_plot, outpath+"allhours/control/control_{date}_{fhour}_sigheight.png".format(date=date, fhour=fhour3))
            fig(coupled_plot, outpath+"allhours/onewaycpl/onewaycpl_{date}_{fhour}_sigheight.png".format(date=date, fhour=fhour3))
            fig(coupled_control_plot, outpath+"allhours/difference/onewaycpl_control_{date}_{fhour}_sigheight_difference.png".format(date=date, fhour=fhour3))
        

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', help="path to output directory", default="./")
    ap.add_argument('-d', '--date', help="date of model run start", required=True)
#    ap.add_argument('-i2', '--input2', help="path to satilite input file", required=True)
#    ap.add_argument('-v', '--variable', help="variable graphed", required=True)
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.date, MyArgs.output)
    #gen_figure(MyArgs.input1, MyArgs.input2, MyArgs.output, MyArgs.variable)


