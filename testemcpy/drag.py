import numpy as np
import matplotlib.pyplot as plt
from emcpy.plots.plots import Scatter
from emcpy.plots.plots import LinePlot
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded
import netCDF4 as nc
import argparse
import os

def gen_figure(inpath, date, outpath):
    def correct(data, variable):
        original=np.array(data.variables[variable][:])
        output=np.full(original.shape, np.nan)
        fill_val=data.variables[variable].__dict__['_FillValue']
        type_val=np.array(data.variables['land'][:])
        for i in range(len(time)):
            for j in range(len(y)):
                for k in range(len(x)):
                    if original[i][j][k] == fill_val:
                        continue
                    elif type_val[i][j][k] == 0:
                        output[i][j][k]=original[i][j][k]
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
        plot.add_colorbar(fontsize=12, extend='neither')
        return plot

    def u10tocd(U10):
        cd=np.full(U10.shape, np.nan)
        for i in range(len(U10)):
            vonkar = 0.4
            nu=1e-6
            grav = 9.81
            alph = np.min((0.028,0.0017 * U10[i] - 0.005))
            
            ust1=0
            ust = U10[i]*np.sqrt(0.001)
            ct=0
            while np.absolute(ust1/ust-1.)>0.001:
                ct=ct+1
                ust1=ust
                z0sm=0.11*nu/ust1
                z0rough=alph*ust1**2/grav
                z0=z0sm+z0rough
                CD = ( vonkar / np.log(10/z0) )**2
                ust = U10[i]*np.sqrt(CD)
                cd[i]=CD

                if ct>20:
                    print('COARE Cd calculation not converging...')

        return cd

 
    def scatterplot(v1, v2, name, xlim, ylim, xlab, ylab):
        x=np.linspace(0,xlim,100)
        y=u10tocd(x)
        lp=LinePlot(x,y)
        lp.label=('COARE 3.5')

        sctr1 = Scatter(v1, v2)
        sctr1.density_scatter()
        plot = CreatePlot()
        plot.plot_layers = [sctr1, lp]
        #plot.plot_layers = [lp]
        plot.add_title(label=name)
        plot.add_xlabel(xlabel=xlab)
        plot.add_ylabel(ylabel=ylab)
        plot.add_legend()
        plot.set_ylim([0,ylim])
        plot.set_xlim([0,xlim])
        return plot

    def fig(data, name):
            fig1 = CreateFigure()
            fig1.plot_list = [data]
            fig1.create_figure()
            fig1.save_figure(name)
            fig1.close_figure()


    control_total_wind_list=np.array([])
    control_total_cd_list=np.array([])
    control_total_sfcr_list=np.array([])
    onewaycpl_total_wind_list=np.array([])
    onewaycpl_total_cd_list=np.array([])
    onewaycpl_total_sfcr_list=np.array([])

    hour_list=[3,6,9,12,18,24,36,48,60,72,96,120,144,168,240,288,336,384]
   
    for hour in hour_list:
        if len(str(hour))==1:
            i3="00"+str(hour)
        elif len(str(hour))==2:
            i3="0"+str(hour)
        else:
            i3=str(hour)
        f = "{inpath}/control/control{date}00/COMROOT/control{date}00/gfs.{date}/00/model_data/atmos/history/gfs.t00z.sfcf{hour}.nc".format(inpath=inpath, hour=i3, date=date)
        control_data=nc.Dataset(f)
        latitude=np.array(control_data.variables['lat'][:])
        longitude=np.array(control_data.variables['lon'][:])
        control_ffmm=np.array(control_data.variables['ffmm'][:])
        control_ugrd10m=np.array(control_data.variables['ugrd10m'][:])
        control_vgrd10m=np.array(control_data.variables['vgrd10m'][:])
        x=np.array(control_data.variables['grid_xt'][:])
        y=np.array(control_data.variables['grid_yt'][:])
        time=np.array(control_data.variables['time'][:])
        control_fill_val=control_data.variables['ffmm'].__dict__['_FillValue']
        control_type_val=np.array(control_data.variables['land'][:])
        control_sfcr=np.array(control_data.variables['sfcr'][:])
        
        control_cd=0.16/(correct(control_data, 'ffmm')**2)
        control_sfcr=correct(control_data, 'sfcr')
        control_u=correct(control_data, 'ugrd10m')
        control_v=correct(control_data, 'vgrd10m')
        control_wind_speed=np.sqrt(control_u**2+control_v**2)
        control_wind_list=np.ravel(control_wind_speed)
        control_cd_list= np.ravel(control_cd)
        control_sfcr_list=np.ravel(control_sfcr)
        control_total_wind_list=np.append(control_total_wind_list, control_wind_list)
        control_total_cd_list=np.append(control_total_cd_list, control_cd_list)
        control_total_sfcr_list=np.append(control_total_sfcr_list, control_sfcr_list)
        
        label="{model} {date} forecast hour {fhour}".format(model='Two-way coupled', date=date, fhour=i3)
        drag_vs_wind=scatterplot(control_wind_list, control_cd_list, label, 35, .003, 'Wind(m/s)', 'Drag coefficient')
        drag_vs_sfcr=scatterplot(control_sfcr_list, control_cd_list, label, .004, .003, 'Surface roughness', 'Drag coefficient')
        sfcr_vs_wind=scatterplot(control_wind_list, control_sfcr_list, label, 35, .004, 'Wind(m/s)', 'Surface roughness')
        drag_map=mapplot(control_cd[0].transpose(), 0, .003, 'plasma', label+' drag coefficient')
        wind_map=mapplot(control_wind_speed[0].transpose(), 0, 35, 'turbo', label+' wind speed')
        sfcr_map=mapplot(control_sfcr[0].transpose(), 0, .004, 'turbo', label+' surface roughness')
        
        model='control'
        fig(drag_map, outpath+'{model}_{date}_{fhour}_drag_coefficient_map.png'.format(model=model, date=date, fhour=i3))
        fig(wind_map, outpath+'{model}_{date}_{fhour}_wind_map.png'.format(model=model, date=date, fhour=i3))
        fig(sfcr_map, outpath+'{model}_{date}_{fhour}_sfcr_map.png'.format(model=model, date=date, fhour=i3))
        fig(drag_vs_wind, outpath+'{model}_{date}_{fhour}_drag_vs_wind.png'.format(model=model, date=date, fhour=i3))
        fig(drag_vs_sfcr, outpath+'{model}_{date}_{fhour}_drag_vs_sfcr.png'.format(model=model, date=date, fhour=i3))
        fig(sfcr_vs_wind, outpath+'{model}_{date}_{fhour}_sfcr_vs_wind.png'.format(model=model, date=date, fhour=i3))
        
        f = "{inpath}/{dtype}/{dtype}{date}00/COMROOT/{dtype}{date}00/gfs.{date}/00/model_data/atmos/history/gfs.t00z.sfcf{hour}.nc".format(inpath=inpath, dtype='onewaycpl', hour=i3, date=date)
        onewaycpl_data=nc.Dataset(f)
        onewaycpl_ffmm=np.array(onewaycpl_data.variables['ffmm'][:])
        ugrd10m=np.array(onewaycpl_data.variables['ugrd10m'][:])
        vgrd10m=np.array(onewaycpl_data.variables['vgrd10m'][:])
        x=np.array(onewaycpl_data.variables['grid_xt'][:])
        y=np.array(onewaycpl_data.variables['grid_yt'][:])
        latitude=np.array(onewaycpl_data.variables['lat'][:])
        longitude=np.array(onewaycpl_data.variables['lon'][:])
        time=np.array(onewaycpl_data.variables['time'][:])
        print(np.shape(time))
        fill_val=onewaycpl_data.variables['ffmm'].__dict__['_FillValue']
        type_val=np.array(onewaycpl_data.variables['land'][:])
        sfcr=np.array(onewaycpl_data.variables['sfcr'][:])
        
        onewaycpl_cd=0.16/(correct(onewaycpl_data, 'ffmm')**2)
        onewaycpl_sfcr=correct(onewaycpl_data, 'sfcr')
        onewaycpl_u=correct(onewaycpl_data, 'ugrd10m')
        onewaycpl_v=correct(onewaycpl_data, 'vgrd10m')
        onewaycpl_wind_speed=np.sqrt(onewaycpl_u**2+onewaycpl_v**2)
        onewaycpl_wind_list=np.ravel(onewaycpl_wind_speed)
        onewaycpl_cd_list= np.ravel(onewaycpl_cd)
        onewaycpl_sfcr_list=np.ravel(onewaycpl_sfcr)
        onewaycpl_total_wind_list=np.append(onewaycpl_total_wind_list, onewaycpl_wind_list)
        onewaycpl_total_cd_list=np.append(onewaycpl_total_cd_list, onewaycpl_cd_list)
        onewaycpl_total_sfcr_list=np.append(onewaycpl_total_sfcr_list, onewaycpl_sfcr_list)
        
        model='onewaycpl'
        label="One-way coupled {date} forecast hour {fhour}".format(model=model, date=date, fhour=i3)
        drag_vs_wind=scatterplot(onewaycpl_wind_list, onewaycpl_cd_list, label, 35, .003, 'Wind(m/s)', 'Drag coefficient')
        drag_vs_sfcr=scatterplot(onewaycpl_sfcr_list, onewaycpl_cd_list, label, .004, .003, 'Surface roughness', 'Drag coefficient')
        sfcr_vs_wind=scatterplot(onewaycpl_wind_list, onewaycpl_sfcr_list, label, 35, .004, 'Wind(m/s)', 'Surface roughness')
        drag_map=mapplot(onewaycpl_cd[0].transpose(), 0, .003, 'plasma', label+' drag coefficient')
        wind_map=mapplot(onewaycpl_wind_speed[0].transpose(), 0, 35, 'turbo', label+' wind speed')
        sfcr_map=mapplot(onewaycpl_sfcr[0].transpose(), 0, .004, 'turbo', label+' surface roughness')

        fig(drag_map, outpath+'{model}_{date}_{fhour}_drag_coefficient_map.png'.format(model=model, date=date, fhour=i3))
        fig(wind_map, outpath+'{model}_{date}_{fhour}_wind_map.png'.format(model=model, date=date, fhour=i3))
        fig(sfcr_map, outpath+'{model}_{date}_{fhour}_sfcr_map.png'.format(model=model, date=date, fhour=i3))
        fig(drag_vs_wind, outpath+'{model}_{date}_{fhour}_drag_vs_wind.png'.format(model=model, date=date, fhour=i3))
        fig(drag_vs_sfcr, outpath+'{model}_{date}_{fhour}_drag_vs_sfcr.png'.format(model=model, date=date, fhour=i3))
        fig(sfcr_vs_wind, outpath+'{model}_{date}_{fhour}_sfcr_vs_wind.png'.format(model=model, date=date, fhour=i3))

        
        cd_diff=onewaycpl_cd-control_cd
        sfcr_diff=onewaycpl_sfcr-control_sfcr
        cd_diffmax=np.nanmax(np.absolute(cd_diff))
        sfcr_diffmax=np.nanmax(np.absolute(sfcr_diff))
        cd_diff_plot=mapplot(cd_diff[0].transpose(), -1*cd_diffmax, cd_diffmax, 'seismic', 'One-way coupled - two-way coupled {date} forecast hour {fhour} difference in drag coefficient'.format(date=date, fhour=i3))
        sfcr_diff_plot=mapplot(sfcr_diff[0].transpose(), -1*sfcr_diffmax, sfcr_diffmax, 'seismic', 'One-way coupled - two-way coupled {date} forecast hour {fhour} difference in surface roughness'.format(date=date, fhour=i3))
        fig(cd_diff_plot, outpath+'{date}_{fhour}_drag_coefficient_difference_map.png'.format(date=date, fhour=i3))
        fig(sfcr_diff_plot, outpath+'{date}_{fhour}_sfcr_difference_map.png'.format(date=date, fhour=i3))
            
                         
    control_total_drag_vs_wind=scatterplot(control_total_wind_list, control_total_cd_list, '{model} drag coefficient vs wind'.format(model='Two-way coupled'), 35, .003, 'Wind (m/s)', 'Drag coefficient')
    control_total_sfcr_vs_wind=scatterplot(control_total_wind_list, control_total_sfcr_list, '{model} sfcr vs wind'.format(model='Two-way coupled'), 35, .004, 'Wind (m/s)', 'Surface roughness')

    onewaycpl_total_drag_vs_wind=scatterplot(onewaycpl_total_wind_list, onewaycpl_total_cd_list, '{model} drag coefficient vs wind'.format(model='One-way coupled'), 35, .003, 'Wind (m/s)', 'Drag coefficient')
    onewaycpl_total_sfcr_vs_wind=scatterplot(onewaycpl_total_wind_list, onewaycpl_total_sfcr_list, '{model} surface roughness vs wind'.format(model='One-way coupled'), 35, .004, 'Wind (m/s)', 'Surface roughness')


    fig(control_total_drag_vs_wind, outpath+'{model}_{date}_overall_drag_vs_wind.png'.format(model='control', date=date))
    fig(onewaycpl_total_drag_vs_wind, outpath+'onewaycpl_{date}_overall_drag_vs_wind.png'.format(date=date))
    fig(control_total_sfcr_vs_wind, outpath+'{model}_{date}_overall_sfcr_vs_wind.png'.format(model='control', date=date))
    fig(onewaycpl_total_sfcr_vs_wind, outpath+'{model}_{date}_overall_sfcr_vs_wind.png'.format(model='onewaycpl', date=date))



if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', help="start of path to input files, does not include model type", required=True)
    #ap.add_argument('-m', '--model', help="model type", required=True)
    ap.add_argument('-d', '--date', help="initial date of run", required=True)
    ap.add_argument('-o', '--output', help="outpath", required=True)
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.date, MyArgs.output)



