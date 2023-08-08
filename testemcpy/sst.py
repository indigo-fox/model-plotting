import numpy as np
import matplotlib.pyplot as plt
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded
from emcpy.plots.plots import LinePlot
import netCDF4 as nc
import argparse
import os


def gen_figure(date):

    month=str(date)[5]
    print(month)
    #Replace fill values with nan, returns corrected values, latitude, longitude
    def correct(dataset, variable, lats, lons):
        fill_val = dataset.variables[variable].__dict__['_FillValue']
        units = dataset.variables[variable].__dict__['units']
        name = dataset.variables[variable].__dict__['long_name']
        land = np.array(dataset.variables['LAND_surface'][:])

        original = np.array(dataset.variables[variable][:])
        corrected = np.full(original[0].shape, np.nan)
        n=0
        for i in range(len(lats)):
            for j in range(len(lons)):
                if original[0][i][j] == fill_val:
                    continue
                elif land[0][i][j] == 0:
                    corrected[i][j] = original[0][i][j]
                else:
                    continue
        return corrected

    def model_avg(modtype):
        model_values=np.full((128, 721, 1440), np.nan)
        for i in range(128):
            hour= 3*i
            if len(str(hour))==1:
                i3="00"+str(hour)
            elif len(str(hour))==2:
                i3="0"+str(hour)
            else:
                i3=str(hour)
            f = "/work2/noaa/marine/ifox/{modtype}/{modtype}{date}00/COMROOT/{modtype}{date}00/gfs.{date}/00/model_data/atmos/master/gfs.t00z.sfluxgrbf{fhour}.grib2.nc".format(fhour=i3, modtype=modtype, date=date)
            dataset=nc.Dataset(f)
            modlats = np.array(dataset.variables['latitude'][:])
            modlons = np.array(dataset.variables['longitude'][:])
            
            corrected=correct(dataset, 'FDNSSTMP_surface', modlats, modlons)
            model_values[i]=corrected
        model_daily=np.full((16, 721, 1440), np.nan)
        for i in range(16):
            model_daily[i]=np.nanmean(model_values[8*i:8*i+8], axis=0)
        return model_daily, modlats, modlons

    onewaycpl_daily, lats, lons=model_avg('onewaycpl')
    control_daily, lats, lons=model_avg('control')
    analysis_daily=np.full((16, 721, 1440), np.nan)
    for d in range(16):
        day=d+13
        f = '/work2/noaa/marine/ifox/Analysis/OSTIA/sst_OSTIA.20200{month}{day}12.0p25.nc'.format(month=month, day=day)
        dataset=nc.Dataset(f)
        alats = np.array(dataset.variables['lat'][:])
        alons = np.array(dataset.variables['lon'][:])
        fill_val = dataset.variables['analysed_sst'].__dict__['_FillValue']
        analysis_sst = np.array(dataset.variables['analysed_sst'][:])
        for i in range(len(alats)):
            for j in range(len(alons)):
                if analysis_sst[0][i][j] == fill_val:
                    analysis_sst[0][i][j] = np.nan
                else:
                    continue
        analysis_daily[d]=analysis_sst


    onewaycpl_analysis_diff=onewaycpl_daily-analysis_daily
    control_analysis_diff=control_daily-analysis_daily
    onewaycpl_control_diff=onewaycpl_daily-control_daily

    control, lats, lons=model_avg('control')

    X, Y = np.meshgrid(lats, lons)
    def map(data, minimum, maximum, cmap, name):
        # Create gridded map object
        gridded = MapGridded(X, Y, np.transpose(data))
        gridded.cmap = cmap
        gridded.vmin = minimum
        gridded.vmax = maximum

        # Create plot object and add features
        plot1 = CreatePlot()
        plot1.plot_layers = [gridded]
        plot1.projection = 'plcarr'
        plot1.domain = 'global'
        plot1.add_map_features(['coastline'])
        #plot1.add_xlabel(xlabel='longitude')
        #plot1.add_ylabel(ylabel='latitude')
        plot1.add_title(label=name, loc='center', fontsize=14)
        plot1.add_grid()
        plot1.add_colorbar(label='K', fontsize=12, extend='neither')
        return plot1

    def fig(plot, name):
        fig1 = CreateFigure()
        fig1.plot_list = [plot]
        fig1.create_figure()
        fig1.save_figure(name)
        fig1.close_figure()

    print(onewaycpl_analysis_diff.shape)
    
    onewaycpl_bias=np.nanmean(onewaycpl_analysis_diff, axis=(1,2))
    control_bias=np.nanmean(control_analysis_diff, axis=(1,2))
    print(len(onewaycpl_bias), len(control_bias))
    onewaycpl_rmse=np.sqrt(np.nanmean(onewaycpl_analysis_diff**2, axis=(1,2)))
    control_rmse=np.sqrt(np.nanmean(control_analysis_diff**2, axis=(1,2)))
    print(len(onewaycpl_rmse), len(control_rmse))
    
    x=np.arange(1, 17, 1)
    lp1=LinePlot(x, onewaycpl_rmse)
    lp1.label=('One-way coupled RMSE')
    lp2=LinePlot(x, control_rmse)
    lp2.label=('Two-way coupled RMSE')
    lp2.color = 'tab:red'
    rmse_plot=CreatePlot()
    rmse_plot.plot_layers=[lp1, lp2]
    rmse_plot.add_title('{} RMSE of SST'.format(date), fontsize=14)
    rmse_plot.add_xlabel('Forecast day')
    rmse_plot.add_ylabel('RMSE (K)')
    rmse_plot.add_legend(loc='upper right')
    fig(rmse_plot, 'images/sst/{}_sst_rmse.png'.format(date))
    
    lp3=LinePlot(x, onewaycpl_bias)
    lp3.label=('One-way coupled bias')
    lp4=LinePlot(x, control_bias)
    lp4.label=('Two-way coupled bias')
    lp4.color = 'tab:red'
    bias_plot=CreatePlot()
    bias_plot.plot_layers=[lp3, lp4]
    bias_plot.add_title('{} SST bias'.format(date), fontsize=14)
    bias_plot.add_xlabel('Forecast day')
    bias_plot.add_ylabel('bias (K)')
    bias_plot.add_legend(loc='upper right')
    fig(bias_plot, 'images/sst/{}_sst_bias.png'.format(date))
    
    
    datamin=np.nanmin((onewaycpl_daily, control_daily, analysis_daily))
    datamax=np.nanmax((onewaycpl_daily, control_daily, analysis_daily))
    adiffmax=np.nanmax(np.absolute((onewaycpl_analysis_diff, control_analysis_diff)))
    moddiffmax=np.nanmax(np.absolute(onewaycpl_control_diff))
    print(adiffmax, moddiffmax)
    print(np.nanpercentile((np.absolute(onewaycpl_analysis_diff), np.absolute(control_analysis_diff)), 99))
    print(np.nanpercentile(np.absolute(onewaycpl_control_diff), 99))
    
    
    for i in range(16):
        if len(str(i+1))==1:
            d='0'+str(i+1)
        else:
            d=str(i+1)
        controlplot=map(control_daily[i], 270, 305, 'turbo', 'Two-way coupled forecast day {day} average SST'.format(day=d))
        onewaycplplot=map(onewaycpl_daily[i], 270, 305, 'turbo', 'One-way coupled forecast day {day} average SST'.format(day=d))
        analysisplot=map(analysis_daily[i], 270, 305, 'turbo', 'Analysis forecast day {day} average SST'.format(day=d))
        fig(controlplot, 'images/sst/control_{date}_{day}_SST.png'.format(date=date, day=d))
        fig(onewaycplplot, 'images/sst/onewaycpl_{date}_{day}_SST.png'.format(date=date, day=d))
        fig(analysisplot, 'images/sst/analysis_{date}_{day}_SST.png'.format(date=date, day=d))


        for j in [6]:
            onewaycpl_analysis_plot=map(onewaycpl_analysis_diff[i], -1*j, j, 'seismic', 'One-way coupled - analysis forecast day {day} average SST'.format(day=d))
            control_analysis_plot=map(control_analysis_diff[i], -1*j, j, 'seismic', 'Two-way coupled - analysis forecast day {day} average SST'.format(day=d))
            onewaycpl_control_plot_samescale=map(onewaycpl_control_diff[i], -1*j, j, 'seismic', 'One-way coupled - Two-way coupled forecast day {day} average SST'.format(day=d))
            fig(onewaycpl_analysis_plot, 'images/sst/onewaycpl_analysis_{date}_{day}_SST_difference_scale{j}.png'.format(date=date, day=d, j=j))
            fig(control_analysis_plot, 'images/sst/control_analysis_{date}_{day}_SST_difference_scale{j}.png'.format(date=date, day=d, j=j))
            fig(onewaycpl_control_plot_samescale, 'images/sst/onewaycpl_control_{date}_{day}_SST_difference_samescale.png'.format(date=date, day=d))


        for k in [2]:
            onewaycpl_control_plot=map(onewaycpl_control_diff[i], -1*k, k, 'seismic', 'One-way coupled - two-way coupled forecast day {day} average SST'.format(day=d))
            fig(onewaycpl_control_plot, 'images/sst/onewaycpl_control_{date}_{day}_SST_difference_scale{k}.png'.format(date=date, day=d, k=k))
    


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    #ap.add_argument('-o', '--output', help="path to output directory", default="./")
    #ap.add_argument('-i1', '--input1', help="path to satilite input file", required=True)
    #ap.add_argument('-i2', '--input2', help="path to satilite input file", required=True)
    ap.add_argument('-d', '--date', help="date of start of run", required=True)
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.date)


