'''
AUTHOR: Kenneth Distefano, June 2026

REPOSITORY: 

PURPOSE: 
create two xth versus m heatmaps to show 
    1) which params tend towards fixation to understand an eradication mechanism. 
    2) segration index
Operations include: 
    - parsing all data files within their respective directories;
    - determing how many runs result in S fixation, coexistence or 
        R fixation;
    - plot results.

NOTE to self:
    ~ need to create general heatmap
    ~ clean up commented out portions
    ~ need to update color bar for segregation index
    ~ save either aveComp, segIndex, or both heat maps


MIT LICENSE:
Copyright <2026> <Kenneth Distefano>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


'''
################################################ imports
import matplotlib.pyplot as plt
import numpy as np
import glob
import argparse
import matplotlib as mpl                # for sequential color map
import re                               # to parse each lattice site
from matplotlib.colors import LinearSegmentedColormap   # create cmap for x plot type
import scipy.io

################################################ definitions
def get_RGB_from_lat_config(infile, comp_thresh):
    '''
    PURPOSE: create a single RGB value that is averaged over every Site
    
    RETURNS: 
        - 1d array with 3 elements: [r, g, b]
        - segregation index of type float
    '''

    # local vars
    rgb = np.zeros(3)
    red = np.array([1.0, 0.0, 0.0])
    black = np.array([0.0, 0.0, 0.0])
    blue = np.array([0.0, 0.0, 1.0])
    seg_index_sum = 0

    # open file, get all data from file, close file
    with open(infile, 'r') as inf:
        rawData = inf.readlines()
    
    # organize data within d=3 array d=(MCS, Lrow, Lcol)
    for line in rawData:
        line = line.strip().split()
        
        # check if new MCS is beginning
        if line[0] == '#':
            mcs = int(line[2])
            row = 0
        
        # get populations of each species and carrying capacity
        else:
            # loop through row of Sites
            for col, Site in enumerate(line):
                # grab populations of each species and carrying capacity with regexes
                # returns list of [N_A, N_B, ... , N_i, K]
                speciesPop = re.findall(r'-?\d+', Site)
                speciesPop = np.array([int(i) for i in speciesPop])

                # save product of R and S composition for segration index later
                seg_index_sum+= (speciesPop[0]*speciesPop[1])/sum(speciesPop[:-1])**2

                # determine which "last_lat_*.dat" color scheme
                if PLOTTYPE == 'Ni':
                    # determine if site is coexisting
                    if (speciesPop > 0).all():
                        rgb += black
                    
                    # if not, check if there's no S --> blue
                    elif speciesPop[0] == 0:
                        rgb += blue

                    # must be zero R --> red
                    else:
                        rgb += red

                # plot type must be "x" --> compute composition for each site
                else:
                    rgb += np.array([speciesPop[0]/sum(speciesPop[:-1]),
                                    0.0,
                                    speciesPop[1]/sum(speciesPop[:-1])])
                    
            
            # update row 
            row += 1

    # compute segregation index
    seg_index = 1 - (1/(L*L))*seg_index_sum/(comp_thresh*(1-comp_thresh))
    

    return rgb/(L*L), seg_index


def indicate_row(y, xi, xf, dx, dy, ax=None, **kwargs):
    '''
    pass desired bottleneck as a row number and return a border around the entire row

    NOTE: Dec 2024, indicate_row() is depreciated --> replaced by horizontal line
    '''
    # create rectangle
    # rect = plt.Rectangle((xi-(dx/2), y-(dy/2)), xf, dy,
    #                      fill=False, **kwargs)
    rect = plt.Rectangle((xi, y-(dy/2)), xf+(dx/2), dy,
                         fill=False, **kwargs)
    ax = ax or plt.gca()
    ax.add_patch(rect)

    

    return rect

def save_avgSegIndex(nestedList):
    '''
    PURPOSE:
        - compute average and standard error of segregation index across realizations
        - save data as *.mat file for exterior MATLAB plotting script
    '''
    # local vars
    

    # loop through composition thresholds
    for row, xth in enumerate(XTH):

        # store data
        avgSegIndexWithStderr = np.zeros((len(M), 3))
        
        # loop through migrations
        for col, m in enumerate(M):

            # compute average across realizations
            avgSegIndexWithStderr[col]=compute_avg_stderr(nestedList[row][col],float(m))

        # save
        # '''testing'''
        # print(f'within save_avgSegIndex():\txth={xth}\n{avgSegIndexWithStderr}')
        save_avgSegIndex_toMatFile(avgSegIndexWithStderr, xth)

    return

def save_avgSegIndex_toMatFile(data2save, xth):
    '''
    PURPOSE:

    passed arguments:
        data2save: len(M)x3 np.array
            - 1st column contains migration rates
            - 2nd col is averaged segration index across realizations
            - 3rd col is standard error of segration index across realizations
    '''

    # create output name
    outputName = f'xth{xth}/segIndex_versus_D_last_L{L}_mcs{MCS}_a0.25_s0.1_K{K[0]}-{K[-1]}_xth{xth}_nu{NU}_delta{DELTA}.mat'

    # save 
    scipy.io.savemat(outputName,
                     {"coexistSites":-1*np.ones((len(M), 3)),
                      "allSites":data2save})

    return

def compute_avg_stderr(dataList, m):
    '''
    PURPOSE:

    passed arguments:
        - dataList: list of data of length number of realizations for a particular 
            composition threshold and migration rate

        - m: migration rate. Needs to be of type float.
    
    '''
    # local vars
    avg = np.average(dataList)
    stderr = np.std(dataList)/np.sqrt(len(dataList))

    # '''testing'''
    # print(f'\twithin compute_avg_stderr():\n\t{dataList}\n\tm= {m}\tavg= {avg}\tstderr= {stderr}')
    

    return np.array([m, avg, stderr])

################################################ main
# create parser
parser = argparse.ArgumentParser(description="create two xth versus m heatmaps to show (1) which params tend towards fixation to understand an eradication mechanism and (2) the computed segration index. Operations include: (i) parsing all data files within their respective directories; (ii) determing how many runs result in S fixation, coexistence or R fixation; (iii) plot results.")

parser.add_argument('-xth', "--threshold", type=float, nargs='*', 
                    default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                    help='(optional flag) x threshold. default is [0.1, ..., 0.9].')
parser.add_argument('-m', '--migs', nargs='*', type=str,
                    default=['0', '0.00001', '0.0000316228', '0.0001',
                             '0.000316228', '0.001', '0.00316228', '0.01',
                             '0.0316228', '0.1', '0.316228', '1', '3.16228', '10'],
                    help='(optional flag) migration values. default is [0, 10^-5, 10^-4.5, ..., 10^-1] of type string to avoid issues parsing data files.')
parser.add_argument('-save', '--saveFigs', default=False, action='store_true',
                    help='(optional flag) include this flag to save figure(s). default is to show to screen.')
# parser.add_argument('-save', '--saveFigs', choices=['aveComp', 'segIndex', 'both'],
#                     help='(optional flag) save aveComp, segIndex, or both heat maps. default is to show to screen.')
parser.add_argument('-c', '--colorbar', default=False, action='store_true',
                    help='(optional flag) include the color bar with heat map')
parser.add_argument("-nol", "--labels", default=True, action='store_false',
                    help='(optional flag) remove x and y labels if provided')
parser.add_argument("-notl", "--noTickLabels", default=False, action='store_true',
                    help='(optional flag) remove x and y tick labels if provided')
parser.add_argument('-t', '--title', default=False, action='store_true',
                    help='(optional flag) include title to heat map. Default=False')
parser.add_argument('-K', '--carryCaps', type=int, nargs='*', default=[80, 1000],
                    help='(optional flag) value for the harsh and mild carrying capacity. Default is [80, 1000]')
parser.add_argument('-L', '--latSideLength', type=int, default=10,
                    help='(optional flag) side length of lattice. default=10')
parser.add_argument("-n", "--nu",
                    help='nu: environmental switching rate')
parser.add_argument("-d", "--delta",
                    help='delta: environmental bias value')
# parser.add_argument("-b", "--bottleNeck", type=float,
#                     help='desired bottle neck ratio')
parser.add_argument('-mcs', '--monteCarloSteps', type=int,
                    help='number of total Monte Carlo Steps, used to correctly identify \'num_den_*.dat\' files.')
parser.add_argument('-ts', "--timeSteps", type=int, nargs='*',
                    help='deisred time step(s)')
parser.add_argument('-df', '--dataFileType', choices=['num', 'last'],
                    help='type of data file to parse: \"num_*.dat\" or \"last_*.dat\"')
parser.add_argument('-p', '--plotType', choices=['x', 'Ni'],
                    help='type of plot that computes x=N_R/N per site or just uses N_i. Intended to be used with data file type \"last_lat_*.dat\"')


# save command line arguments
args = parser.parse_args()
NU = args.nu
DELTA = args.delta
# desiredBotRatio = args.bottleNeck       # desired bottleneck ratio
desiredTimeSteps = args.timeSteps       # desired time step
SAVE = args.saveFigs                    # bool to save/show fig(s)
XTH = args.threshold                    # list of xth values
M = args.migs                           # list of migration values
K = args.carryCaps                      # list of harsh and mild carrying capacities
CBAR = args.colorbar                    # create color bar for heat map
LABELS = args.labels                    # bool to determine if axis labels are used
NOTICKLABELS = args.noTickLabels        # bool to determine if tick labels are used
TITLE = args.title                      # include title for heat map
L = args.latSideLength                  # lattice side length
MCS = args.monteCarloSteps              # number of monte carlo steps
DFTYPE = args.dataFileType              # type of data file to parse
PLOTTYPE = args.plotType                # compute x or use N_i for each site
TEXTSIZE = 25                           # size for text w/n plot
WIDTH = 3                               # size for tick and spine width
LENGTH = 8                              # ticksize length
ALPHA = 0.15                            # fill_between to show confidence interval

# additional variables
markers = ['o', '^', 's', 'X', 'd']
cmap = mpl.cm.get_cmap('plasma')
colors = cmap(np.linspace(0,0.75,len(desiredTimeSteps)))# max at 0.75 to avoid yellow
# nCompThresh = len(XTH)
# nD = len(M)
xthVSm = np.zeros((len(desiredTimeSteps), len(XTH), len(M), 3)) # rgb vals for hmap
avgSegIndexVSm = np.zeros((len(desiredTimeSteps), len(XTH), len(M))) # avg seg index
allSegIndexVSm = [[[] for c in range(len(M))] for r in range(len(XTH))] # all segindx
dirsDict = {}                           # dict to hold all data to help with parsing
botList = []                            # list of bottlenecks K_+/K_-
Dlist = []                              # list of diffusion rates
horizontalD=np.zeros((len(desiredTimeSteps),len(M),2))# 3d arr 1st column is avg and 
                                            # 2nd column is standard deviation
                                            # of the binonal distribution
# allDirs = glob.glob("Nth40/L20_mcs500_a0.25_s0.1_D0*_K80-*/")
# # allDirs = glob.glob("L20_mcs500_a0.25_s0.1_D0*_K80-*/")
# dirs=[]
# for d in allDirs:
#     if 'K80-64000' in d or 'K80-128000' in d:
#         pass
#     else:
#         dirs.append(d)

'''testing'''
print(f'L= {L}, NU= {NU}, DELTA= {DELTA}, timeStep(s)= {desiredTimeSteps}, SAVE= {SAVE}, df type= {DFTYPE}, plot type= {PLOTTYPE}')
print(f'xth= {XTH}')
print(f'm= {M}')


# loop through composition thresholds
for row, xth in enumerate(sorted(XTH, reverse=True)):

    '''testing'''
    print(f'xth= {xth}')
    
    # loop through migration values
    for col, m in enumerate(M):

        # '''testing'''
        # print(f'\tm= {m}')

        # get directory
        dir = f"xth{xth}/L{L}_mcs{MCS}_a0.25_s0.1_D{m}_K{K[0]}-{K[-1]}/"

        # # add comp thresh and migration value to dictionary
        # if xth not in dirsDict:
        #     dirsDict[xth] = {m : [[] for i in range(len(desiredTimeSteps))]}
        # else:
        #     dirsDict[xth][m] = [[] for i in range(len(desiredTimeSteps))]

        # get data files
        if DFTYPE == 'num':
            dfnames=dir+f'num_den_L{L}_xth{xth}_D{m}_K{K[0]}-{K[-1]}_nu{NU}_delta{DELTA}_*.dat'
            
        else:
            dfnames=dir+f'last_lat_config_L{L}_mcs{MCS}_a0.25_s0.1_D{m}_K{K[0]}-{K[-1]}_*th{xth}_nu{NU}_delta{DELTA}_id*.dat'
            
        
        # glob all data files
        dfs=glob.glob(dfnames)
        
        # make sure data files were globbed
        if not dfs:
            print(f'\n\tWARNING! could not find data files... dfs= {dfs}\n')
            print(f'\t{dfnames}\n')
            exit(1)
        elif len(dfs) is not 50:
            print(f'NOTE: len(dfs)= {len(dfs)} =/= 50')

        # # an rgb arry for every desired time step to increment later
        # rgb = np.zeros((len(desiredTimeSteps), 3))
        
        # loop through each data file
        for df in dfs:

            # parsing "num_*.dat" data files
            if DFTYPE == 'num':
                # get data; returns 2d arr w/ dims=(MCS,2) where cols are S R pops
                data = np.loadtxt(df)


                # loop through desired time slices
                for num, ts in enumerate(desiredTimeSteps):

                    # compute the relative fraction of particles
                    y = int((data[ts,1] - data[ts,0])/np.sum(data[ts]))

                    # '''testing'''
                    # if m == 0 and y is not 0:
                    #     print(f'\ty= {y}, [S,R]= {data[ts]}')

                    # increment rgb value to average later
                    rgb = y*np.array([y-1, 0, y+1])/2
                    xthVSm[num,row,col] += rgb

                    # '''testing'''
                    # print(f'\t\ty= {y}\trgb= {rgb}\nxthVSm[]=\n{xthVSm}\n\n')

            # parsing "last_lat_*.dat" data files
            else:
                # create empty 3d array to hold data from file
                lattice = np.zeros((2, L, L))

                # get rgb data and computed segration index
                rgb, segIndex = get_RGB_from_lat_config(df, xth)

                # save in 4d array
                xthVSm[0,row, col] += rgb

                # save in 3d array
                avgSegIndexVSm[0,row,col] += segIndex

                # save segregation index into list to save in *.mat file later
                allSegIndexVSm[row][col].append(segIndex)



# average rgb value (divide by number of realizations)
xthVSm /= len(dfs)

# average segregation index (divide by number of realizations)
avgSegIndexVSm /= len(dfs)
'''testing'''
print(f'\nsegration index')
print(avgSegIndexVSm)

# save *.mat files for MATLAB plotting script
save_avgSegIndex(allSegIndexVSm)


# # iterate through dictionary of dictionaries
# for row, (bot, DandRuns) in enumerate(sorted(dirsDict.items(), reverse=True)):
    
#     for col, (D, run) in enumerate(sorted(dirsDict[bot].items())):

#         # save bottleneck and diffusion values
#         if bot not in botList:
#             botList.insert(0, bot)
        
#         if D not in Dlist:
#             Dlist.append(D)
        
#         # loop through each desired time step
#         for i, ts in enumerate(desiredTimeSteps):
#             # average runs and save in single array
#             avgBotvsD[i, row, col] = np.average(run[i])

#             # make sure there's 200 realizations per square
#             if len(run[i]) != 200:    
#                 print(f'number of runs= {len(run[i])}\tbottleneck= {bot}\tm= {D}')

#             # for particular bottleneck, compute average with standard deviation
#             if bot == desiredBotRatio:
#                 horizontalD[i,col,0] = np.average(run[i])
                
#                 # computing error --> see confidence intervals, wald method w/ z=1
#                 horizontalD[i,col,1]=np.sqrt(horizontalD[i,col,0]*(1-horizontalD[i,col,0])/len(run[i]))



# '''testing'''
# print(botList, Dlist)
# print(avgBotvsD)
# print('desired bot ratio= {}'.format(desiredBotRatio), horizontalD, sep='\n')

# determine where the x and y tick labels will be
centers = [0, len(M), 0, len(XTH)]
dx, = np.diff(centers[:2]) / (len(M)-1)
dy, = np.diff(centers[2:]) / (len(XTH)-1)
myExtent = [centers[0]-dx/2, centers[1]+dx/2, centers[2]-dy/2, centers[3]+dy/2] #lrbt
xTickLoc = np.arange(centers[0], centers[1]+dx, dx)
yTickLoc = np.arange(centers[2], centers[3]+dy, dy)


# # compute theoretical bottle neck w.r.t. migration
# theoBotNeck = compute_theoBotNeck(Dlist[1:])
# newTheoBotNeck = 2*1.1*np.log2(float(NU)/(12.5*80*np.array([10**-5.25]+Dlist[1:])))

'''testing'''
print(f'centers= {centers}, dx= {dx}, dy= {dy}')
print(f'myExtent= {myExtent}')
print(f'xTickLoc= {xTickLoc}')
print(f'yTickLoc= {yTickLoc}')

# print()
# print(f'mRates= {Dlist[1:]}')

# print()
# print(f'transformed theoBotNeck= {theoBotNeck}')
# print(f'transformed newTheoBotNeck= {newTheoBotNeck}')
# print()

for i, ts in enumerate(desiredTimeSteps):
    # create heat map
    fig, ax = plt.subplots(figsize=(9,6))
    fig_segIndx, ax_segIndx = plt.subplots(figsize=(9,6))

    # plot segration index
    img_segIndx = ax_segIndx.imshow(avgSegIndexVSm[i], vmin=0, vmax=1,
                                    extent=myExtent, aspect='auto', cmap='gray_r')

    # plot 2d array
    if PLOTTYPE == 'x':
        blackBlue = [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
        mycmap_name = 'black_blue'
        xcmap = LinearSegmentedColormap.from_list(mycmap_name, blackBlue)
        img = ax.imshow(xthVSm[i], vmin=0, vmax=1, extent=myExtent, aspect='auto',
                        cmap=xcmap)
    else:
        img = ax.imshow(xthVSm[i], vmin=0, vmax=1, extent=myExtent, aspect='auto')

    # # plot theorectical prediction
    # ax.plot(xTickLoc[1:], theoBotNeck, color='green', linewidth=WIDTH+2)
    # ax.plot(np.insert(xTickLoc[1:],0,dx/2), newTheoBotNeck,
    #         color='gold', linewidth=WIDTH+2)

    # loop through both heat maps
    for a in [ax, ax_segIndx]:
        
        # plot horizontal line denoting a broken axis
        a.vlines(xTickLoc[0]+(dx/2), ymin=myExtent[2], ymax=myExtent[3], 
                 colors='white',
                 linewidth=WIDTH, linestyles='dashed')
    

        # format plot
        if LABELS:
            a.set_xlabel('$m$', fontsize=TEXTSIZE)
            a.set_ylabel('$x_{th}$', fontsize=TEXTSIZE)

        if TITLE:
            a.set_title(f'($\\nu$= {NU}, $\\delta$= {DELTA}) K$\in[{K[0]}, {K[-1]}]$ t={ts}', fontsize=TEXTSIZE)

        # increase the thickness of figure border
        for spine in ['left', 'right', 'top', 'bottom']:    # incr border
            a.spines[spine].set_linewidth(WIDTH)
    
        a.set_ylim(bottom=myExtent[2], top=myExtent[3])    # avoid shrinking heat map
        a.set_xticks(xTickLoc)
        a.set_yticks(yTickLoc)
        a.tick_params(axis='both', which='major', labelsize=TEXTSIZE-3, width=WIDTH,
                        length=LENGTH)

        if NOTICKLABELS:
            a.set_xticklabels([])
            a.set_yticklabels([])

        else:
            # check if migration values go until 10^1 (const env) or 10^-1 (dyn env)
            if len(M) == 14:
                xTickLabels = [f"$10^{{{i:.0f}}}$" for i in np.arange(-5,1.5,0.5)]
            elif len(M) == 10:
                xTickLabels = [f"$10^{{{i:.0f}}}$" for i in np.arange(-5,-0.5,0.5)]
            else:
                print(f'WARNING!! Double-check desired migration values; Could not create xTickLabels.\nm= {M}')
                exit(1)

            # yTickLabels = XTH

            # # remove unnecessary decimals to only select y tick labels
            # for i,y in enumerate(yTickLabels):
            #     if i>0 and i%2 == 0:
            #         yTickLabels[i] = int(y)

            a.set_xticklabels([f'{float(M[0]):.0f}']+xTickLabels)
            a.set_yticklabels(XTH)

            # set every other label to be not visible
            for xl in a.xaxis.get_ticklabels()[2::2]:
                xl.set_visible(False)
                
            # for yl in ax.yaxis.get_ticklabels()[1::2]:
            #     yl.set_visible(False)


    if CBAR:
        # fixation heat map
        cbar = fig.colorbar(mappable=img, ax=ax, ticks=[0, 1])
        cbar.ax.tick_params(labelsize=TEXTSIZE, width=WIDTH, length=LENGTH)
        cbar.ax.set_title('$x=N_R/N$',fontsize=TEXTSIZE, pad=TEXTSIZE-8)
        cbar.outline.set_linewidth(WIDTH)

        # segration index heat map
        cbar_segIndx = fig_segIndx.colorbar(mappable=img_segIndx, ax=ax_segIndx,
                                            ticks=[0, 1])
        cbar_segIndx.ax_segIndx.tick_params(labelsize=TEXTSIZE, width=WIDTH,
                                            length=LENGTH)
        cbar_segIndx.ax_segIndx.set_title('$\\langle \mathcal{S} \\rangle$',
                                          fontsize=TEXTSIZE,
                                          pad=TEXTSIZE-8)
        cbar_segIndx.outline.set_linewidth(WIDTH)


    fig.tight_layout()
    fig_segIndx.tight_layout()
    
    # detemine output file name depending on parameters
    # for simplicity, this only works if the -notl flag was included
    # doesn't consider of ticklabels or colorbar are desired
    ofnames = []
    if LABELS:
        ofnames.append(f'xth/xthVSm_{DFTYPE}'+(f'-{PLOTTYPE}' if PLOTTYPE=='x' else '')+f'_L{L}_mcs{MCS}_K{K[0]}-{K[-1]}_nu{NU}_delta{DELTA}_ts{ts}.png')
        
        ofnames.append(f'xth/segIndex_{DFTYPE}_L{L}_mcs{MCS}_K{K[0]}-{K[-1]}_nu{NU}_delta{DELTA}_ts{ts}.png')
    else:
        ofnames.append(f'xth/xthVSm_{DFTYPE}'+(f'-{PLOTTYPE}' if PLOTTYPE=='x' else '')+f'_L{L}_mcs{MCS}_K{K[0]}-{K[-1]}_nu{NU}_delta{DELTA}_ts{ts}_noLabels.png')

        ofnames.append(f'xth/segIndex_{DFTYPE}_L{L}_mcs{MCS}_K{K[0]}-{K[-1]}_nu{NU}_delta{DELTA}_ts{ts}_noLabels.png')

    # save
    if SAVE:
        print(f'saved {ofnames[-1]}')
        fig_segIndx.savefig(ofnames[-1])
    else:
        print(f'showing {ofnames}')
        plt.show()


'''testing'''
exit(1)



# for a particular bottleneck ratio, create a figure to show nonmonotinicity
# broken axis was created through this tutorial:
# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html
if desiredBotRatio:
    # new plot visualizing horizontal slice of heat map with broken axis
    # set up figure
    fig, axes = plt.subplots(1, 2, figsize=(10,6), sharey=True,
                             gridspec_kw={'width_ratios': [1,10]})

    # loop through time steps backwards so legend order matches plot order
    for i, ts in enumerate(sorted(desiredTimeSteps, reverse=True), 1):
        # plot errorbars and markers for m=0
        axes[0].errorbar(x=Dlist, y=horizontalD[-i,:,0], yerr=horizontalD[-i,:,1],
                         fmt=markers[-i], capsize=10, markersize=10,color=colors[-i])
        
        # plot with lines and with fill between
        axes[1].plot(Dlist, horizontalD[-i,:,0], label='$t= {}$'.format(ts+1),
                     color=colors[-i])
        axes[1].fill_between(Dlist, horizontalD[-i,:,0]-horizontalD[-i,:,1],
                             horizontalD[-i,:,0]+horizontalD[-i,:,1],
                             color=colors[-i], alpha=ALPHA)
    
    # format both axes
    for ax in axes:
        
        # make all axes and tick labels thicker, more readable
        for spine in ['left', 'right', 'top', 'bottom']:
            ax.spines[spine].set_linewidth(WIDTH)
        ax.tick_params(axis='both', which='major', direction='out',
                       labelsize=TEXTSIZE, width=WIDTH, length=LENGTH)
        ax.grid()
    
    # zoom into axes to see desired data
    axes[0].set_xlim(-1e-6,1e-6)

    # hide the spines between ax and ax2 and ticks for ax2
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)

    # Now, let's turn towards the cut-out slanted lines.
    # We create line objects in axes coordinates, in which (0,0), (0,1),
    # (1,0), and (1,1) are the four corners of the Axes.
    # The slanted lines themselves are markers at those locations, such that the
    # lines keep their angle and position, independent of the Axes size or scale
    # Finally, we need to disable clipping.

    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                  linestyle="none", color='k', mec='k', mew=WIDTH, clip_on=False)
    axes[0].plot([1, 1], [0, 1], transform=axes[0].transAxes, **kwargs)
    axes[1].plot([0, 0], [0, 1], transform=axes[1].transAxes, **kwargs)


    # format plot
    axes[0].set_ylabel('Prob. of R extinction', fontsize=TEXTSIZE)
    axes[0].set_xticks([0])
    axes[0].set_xticklabels(['0'])

    axes[1].tick_params('y', which='major', width=0)    #hide tick between ax and ax2
    axes[1].tick_params('x', which='minor', width=WIDTH-1, length=LENGTH-2)
    axes[1].set_xlabel('$m$', fontsize=TEXTSIZE)
    axes[1].set_xscale('log', basex=10)
    # ax.set_ylim(bottom=0, top=1)
    axes[1].legend(fontsize=TEXTSIZE)

    fig.tight_layout()

    # determine output file name depending on parameters
    # for simplicity, this only works if the -notl flag was included
    # doesn't consider of ticklabels or colorbar are desired
    if LABELS:
        ofname = f'Nth40/K{desiredBotRatio*80}vsD_L{L}_Nth40_mcs{MCS}_nu{NU}_delta{DELTA}_ts{desiredTimeSteps[0]}-{desiredTimeSteps[-1]}.png'
    else:
        ofname = f'Nth40/K{desiredBotRatio*80}vsD_L{L}_Nth40_mcs{MCS}_nu{NU}_delta{DELTA}_ts{desiredTimeSteps[0]}-{desiredTimeSteps[-1]}_noLabels.png'

    # save
    print(ofname)
    if SAVE:
        fig.savefig(ofname)
    else:
        plt.show()
    


