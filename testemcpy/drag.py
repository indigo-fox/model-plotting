import numpy as np
import matplotlib.pyplot as plt
from emcpy.plots.plots import Scatter
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded
import netCDF4 as nc
import argparse
import os

def gen_figure(model, date, outpath):
    def correct(variable):
        output=np.full(variable.shape, np.nan)
        for i in range(len(time)):
            for j in range(len(y)):
                for k in range(len(x)):
                    if variable[i][j][k] == fill_val:
                        continue
                    elif type_val[i][j][k] == 0:
                        output[i][j][k]=variable[i][j][k]
                    else:
                        continue
        return output
    
    def mapplot(data, minimum, maximum, cmap, name):

        # Create gridded map object
        gridded = MapGridded(latitude, longitude, np.transpose(data))
        gridded.cmap = cmap
        gridded.vmin = minimum
        gridded.vmax = maximum

        # Create plot object and add features
        plot = CreatePlot()
        plot.plot_layers = [gridded]
        plot.projection = 'plcarr'
        plot.domain = 'global'
        plot.add_map_features(['coastline'])
        plot.add_xlabel(xlabel='longitude')
        plot.add_ylabel(ylabel='latitude')
        plot.add_title(label=name, loc='center')
        plot.add_grid()
        plot.add_colorbar(label='m/s', fontsize=12, extend='neither')
        return plot

    def scatterplot(v1, v2, name):
        sctr1 = Scatter(v1, v2)
        sctr1.density_scatter()
        plot = CreatePlot()
        plot.plot_layers = [sctr1]
        plot.add_title(label='{}: Drag coefficient vs wind speed over ocean'.format(name))
        plot.add_xlabel(xlabel='Wind (m/s)')
        plot.add_ylabel(ylabel='Drag coefficient')
        plot.add_legend()
        plot.set_ylim([0,0.003])
        plot.set_xlim([0,35])

        return plot

    total_wind_list=np.array([])
    total_cd_list=np.array([])
    hour_list=[3,6,9,12,18,24,36,48,60,72,96,120,144,168,192,240,288,336,384]
    for hour in hour_list:
        if len(str(hour))==1:
            i3="00"+str(hour)
        elif len(str(hour))==2:
            i3="0"+str(hour)
        else:
            i3=str(hour)
        f = "/work2/noaa/marine/ifox/{dtype}/{dtype}{date}00/COMROOT/{dtype}{date}00/gfs.{date}/00/model_data/atmos/history/gfs.t00z.sfcf{hour}.nc".format(dtype=model, hour=i3, date=date)
        data=nc.Dataset(f)
        ffmm=np.array(data.variables['ffmm'][:])
        ugrd10m=np.array(data.variables['ugrd10m'][:])
        vgrd10m=np.array(data.variables['vgrd10m'][:])
        x=np.array(data.variables['grid_xt'][:])
        y=np.array(data.variables['grid_yt'][:])
        latitude=np.array(data.variables['lat'][:])
        longitude=np.array(data.variables['lon'][:])
        time=np.array(data.variables['time'][:])
        fill_val=data.variables['ffmm'].__dict__['_FillValue']
        type_val=np.array(data.variables['land'][:])

        cd=0.16/(correct(ffmm)**2)
        u=correct(ugrd10m)
        v=correct(vgrd10m)
        wind_speed=np.sqrt(u**2+v**2)
        wind_list=np.ravel(wind_speed)
        cd_list= np.ravel(cd)
        total_wind_list=np.append(total_wind_list, wind_list)
        total_cd_list=np.append(total_cd_list, cd_list)
        
        label="{model} {date} forecast hour {fhour}".format(model=model, date=date, fhour=i3)
        drag_vs_wind=scatterplot(wind_list, cd_list, label)
        drag_map=mapplot(cd[0].transpose(), 0, .003, 'plasma', label+' Drag coefficient')
        wind_map=mapplot(wind_speed[0].transpose(), 0, 35, 'turbo', label+' wind speed')

        fig1 = CreateFigure()
        fig1.plot_list = [drag_map]
        fig1.create_figure()
        fig1.save_figure(outpath+'{model}_{date}_{fhour}_drag_coefficient_map.png'.format(model=model, date=date, fhour=i3))
        fig1.close_figure()

        fig2 = CreateFigure()
        fig2.plot_list = [drag_vs_wind]
        fig2.create_figure()
        fig2.save_figure(outpath+'{model}_{date}_{fhour}_drag_vs_wind.png'.format(model=model, date=date, fhour=i3))
        fig2.close_figure()

        fig3 = CreateFigure()
        fig3.plot_list = [wind_map]
        fig3.create_figure()
        fig3.save_figure(outpath+'{model}_{date}_{fhour}_wind_map.png'.format(model=model, date=date, fhour=i3))
        fig3.close_figure()



                         
    total_drag_vs_wind=scatterplot(total_wind_list, total_cd_list, '{model} {date} overall'.format(model=model, date=date))
    fig3 = CreateFigure()
    fig3.plot_list = [total_drag_vs_wind]
    fig3.create_figure()
    fig3.save_figure(outpath+'{model}_{date}_overall_drag_vs_wind.png'.format(model=model, date=date))
    fig3.close_figure()
    

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-m', '--model', help="model type", required=True)
    ap.add_argument('-d', '--date', help="initial date of run", required=True)
    ap.add_argument('-o', '--output', help="outpath", required=True)
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.model, MyArgs.date, MyArgs.output)



