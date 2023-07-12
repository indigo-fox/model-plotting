import numpy as np
import matplotlib.pyplot as plt
from emcpy.plots.plots import Scatter
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded
from emcpy.stats import get_linear_regression
import netCDF4 as nc
import argparse
import os

def gen_figure(model, date, outpath):
    def correct(variable):
        output=np.full(variable.shape, np.nan)
        for i in range(len(time)):
            for j in range(len(y)):
                for k in range(len(x)):
                    if variable[i][j][k] == fill_val or variable[i][j][k] == 0:
                        continue
                    elif type_val[i][j][k] == 0:
                        output[i][j][k]=variable[i][j][k]
                    else:
                        continue
        return output
    
    hour_list=[24]
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
        fricv=np.array(data.variables['fricv'][:])
        ugrd_hyblev1=np.array(data.variables['ugrd_hyblev1'][:])
        vgrd_hyblev1=np.array(data.variables['vgrd_hyblev1'][:])

        
        cd1=0.16/(correct(ffmm)**2)
        u10=correct(ugrd10m)
        v10=correct(vgrd10m)
        wind_speed=np.sqrt(u10**2+v10**2)
        uzr=correct(ugrd_hyblev1)
        vzr=correct(vgrd_hyblev1)
        ustar=correct(fricv)
        wszr=np.sqrt(uzr**2+vzr**2)
        cd2=ustar**2/wszr**2

        wind_list=np.ravel(wind_speed)
        cd1_list=np.ravel(cd1)
        cd2_list=np.ravel(cd2)

        sctr1 = Scatter(cd1_list, cd2_list)
        sctr1.density_scatter()
        plot1 = CreatePlot()
        plot1.plot_layers = [sctr1]
        plot1.add_title(label=("cd1 vs cd2"))
        plot1.add_xlabel(xlabel='cd1')
        plot1.add_ylabel(ylabel='cd2')
        plot1.add_legend()
        plot1.set_ylim([0,0.003])
        plot1.set_xlim([0,0.003])

        sctr3 = Scatter(cd1_list, cd2_list)
        sctr3.density_scatter()
        plot3 = CreatePlot()
        plot3.plot_layers = [sctr3]
        plot3.add_title(label=("cd1 vs cd2"))
        plot3.add_xlabel(xlabel='cd1')
        plot3.add_ylabel(ylabel='cd2')
        plot3.add_legend()
        plot3.set_ylim([0,1])
        plot3.set_xlim([0,0.003])

        sctr2 = Scatter(cd1_list, cd2_list)
        sctr2.density_scatter()
        plot2 = CreatePlot()
        plot2.plot_layers = [sctr2]
        plot2.add_title(label=("cd1 vs cd2"))
        plot2.add_xlabel(xlabel='cd1')
        plot2.add_ylabel(ylabel='cd2')
        plot2.add_legend()

        fig1 = CreateFigure()
        fig1.plot_list = [plot1]
        fig1.create_figure()
        fig1.save_figure(outpath+'0.003_cd1_vs_cd2_'+i3+'.png')
        fig1.close_figure()

        fig2 = CreateFigure()
        fig2.plot_list = [plot2]
        fig2.create_figure()
        fig2.save_figure(outpath+'cd1_vs_cd2_'+i3+'.png')
        fig2.close_figure()

        fig3 = CreateFigure()
        fig3.plot_list = [plot3]
        fig3.create_figure()
        fig3.save_figure(outpath+'1_cd1_vs_cd2_'+i3+'.png')
        fig3.close_figure()

    
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-m', '--model', help="model type", required=True)
    ap.add_argument('-d', '--date', help="initial date of run", required=True)
    ap.add_argument('-o', '--output', help="outpath", required=True)
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.model, MyArgs.date, MyArgs.output)

