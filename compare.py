import numpy as np
import matplotlib.pyplot as plt
from emcpy.plots import CreatePlot, CreateFigure
from emcpy.plots.map_tools import Domain, MapProjection
from emcpy.plots.map_plots import MapGridded
import netCDF4 as nc
import argparse
import os


def gen_figure(inpath1, inpath2, outpath, variable):
    datanc1 = nc.Dataset(inpath1)
    datanc1 = nc.Dataset(inpath2)
    units = datanc1.variables[variable].__dict__['units']
    name = datanc1.variables[variable].__dict__['long_name']


    def correct(dataset):
        lats = np.array(dataset.variables['latitude'][:])
        lons = np.array(dataset.variables['longitude'][:])
        fill_val = dataset.variables[variable].__dict__['_FillValue']
#        units = dataset.variables[variable].__dict__['units']
#        name = dataset.variables[variable].__dict__['long_name']
        original = np.array(dataset.variables[variable][:])
        corrected = np.full(original[0].shape, np.nan)
        for i in range(len(lats)):
            for j in range(len(lons)):
                if original[0][i][j] == fill_val:
                    continue
                else:
                    corrected[i][j] = original[0][i][j]
        return corrected, lats, lons
    data1, lats, lons = correct(nc.Dataset(inpath1))
    data2, lats, lons = correct(nc.Dataset(inpath2))
    datap01 = min(np.nanpercentile(data1, 1), np.nanpercentile(data2, 1))
    datap99 = max(np.nanpercentile(data1, 99), np.nanpercentile(data2, 99))
    
    diff = data1 - data2
    diffp99 = np.nanpercentile(np.absolute(diff), 99)
    
                
    X, Y = np.meshgrid(lats, lons)
    Z = diff
    def map(data, minimum, maximum, cmap):
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
        plot1.add_xlabel(xlabel='longitude')
        plot1.add_ylabel(ylabel='latitude')
        plot1.add_title(label=name, loc='center')
        plot1.add_grid()
        plot1.add_colorbar(label=units,
                           fontsize=12, extend='neither')
        return plot1

    plot1 = map(data1, datap01, datap99, 'viridis')
    plot2 = map(data2, datap01, datap99, 'viridis')
    plot3 = map(diff, -1*diffp99, diffp99, 'bwr')

    # Create figure
    fig1 = CreateFigure()
    fig1.plot_list = [plot1]
    fig1.create_figure()
    fig1.save_figure('{}_cpl.png'.format(outpath))
    fig1.close_figure()

    fig2 = CreateFigure()
    fig2.plot_list = [plot2]
    fig2.create_figure()
    fig2.save_figure('{}_ctrl.png'.format(outpath))
    fig2.close_figure()

    fig3 = CreateFigure()
    fig3.plot_list = [plot3]
    fig3.create_figure()
    fig3.save_figure('{}_difference.png'.format(outpath))
    fig3.close_figure()



    


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', help="path to output directory", default="./")
    ap.add_argument('-i1', '--input1', help="path to satilite input file", required=True)
    ap.add_argument('-i2', '--input2', help="path to satilite input file", required=True)
    ap.add_argument('-v', '--variable', help="variable graphed", required=True)
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input1, MyArgs.input2, MyArgs.output, MyArgs.variable)
