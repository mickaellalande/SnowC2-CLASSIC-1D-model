import os
import sys
import warnings

import cartopy.crs as ccrs
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import yaml

warnings.filterwarnings("ignore")

# first argument is the Processed_Sites directory
inputFiles=sys.argv[1]
print(inputFiles)

# optional second argument is the output path
if len(sys.argv) > 2:
    outputFiles = sys.argv[2] + '/'
else:
    outputFiles = ''

variables = {}
variables['Latent Heat'] = {'model': 'hfls', 'observations': 'LE_F_MDS_filtered'}
variables['Sensible Heat'] = {'model': 'hfss', 'observations': 'H_F_MDS_filtered'}
variables['Gross Primary Productivity'] = {'model': 'gpp', 'observations': 'GPP_DT_VUT_MEAN_filtered'}
variables['Ecosystem Respiration'] = {'model': 'reco', 'observations': 'RECO_DT_VUT_MEAN_filtered'}
#variables['Net Ecosystem Exchange'] = {'model': 'nep', 'observations': 'NEE_VUT_MEAN_filtered'}

# lists and  their order:
landscapedict = {'agriculture' : 0, 'disturbed': 1, 'undisturbed': 2,  'fire': 3, 'harvested': 4 }
biomedict = {'EBF': 0, 'DBF': 1, 'ENF': 2, 'DNF': 3, 'CSH': 4, 'OSH': 5, 'SAV': 6,'GRA': 7, 'WET': 8, 'CRO': 9, 'MF' : 10  }
peatdict = { 'upland': 'o', 'peatland': 's'}

# Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
#colors = sns.color_palette("hls", len(biomedict),desat = 0.5)
colors2 = sns.color_palette(None, len(biomedict))
colors3 = sns.color_palette()

# Get info from site yamls
def yaml_parse(file, att):
    try:
        with open(file, "r") as f:
            d = yaml.safe_load(f)
        return(d[att])
    except:
        d = ""
        return d

# Get site info
def read_sites():
    # import my list of sites.
    sites = []
    for site in os.listdir(inputFiles):
        if os.path.isdir(inputFiles + "/" + site):
            sites.append(site)

    data = []
    for site in sites:
        print("starting " + site)
        ymlpath = inputFiles + "/" + site + "/siteinfo.yaml"
        row = []
        row.append(site)
        if yaml_parse(ymlpath, "name") == "":
            print(f"Could not parse {site} yaml, skipping...")
            continue
        row.append(yaml_parse(ymlpath, "lat"))
        row.append(yaml_parse(ymlpath, "lon"))
        row.append(yaml_parse(ymlpath, "biome"))
        st = yaml_parse(ymlpath, "start")
        en = yaml_parse(ymlpath, "end")
        totyrs = en - st + 1
        row.append(totyrs)
        row.append(yaml_parse(ymlpath, "landscape_state"))
        peatbool = yaml_parse(ymlpath, "peatland")
        if (peatbool == 'true'):
            row.append("peatland")
        else:
            row.append("upland")
        data.append(row)
    df = pd.DataFrame(data, columns=["Site","latitude","longitude","biome","total_years","landscape_state","peatland"])
    print(df)

    return df

# Build biome plot
def biome_plot(df):
    plt.figure(figsize=(4*4, 4*5))

    plot_crs = ccrs.Robinson(central_longitude=0)
    ax = plt.axes(projection=plot_crs)
    ax.set_global()

    ax.coastlines('110m', linewidth=0.8, zorder=2,color='grey')

    # Sites by biome
    for n,b in enumerate(biomedict):
        bio=df[df['biome'].str.match(b)]
        if b == 'SAV':
            bio.plot.scatter(x='longitude', y='latitude',s=4*bio['total_years'],ax=ax,label=b,color='tab:pink',
                    linewidth=1,  marker='o',alpha=0.9,edgecolors='black',
                    transform=ccrs.PlateCarree()
                    )
        elif b == 'WET':
            bio.plot.scatter(x='longitude', y='latitude',s=4*bio['total_years'],ax=ax,label=b,color='tab:olive',
                    linewidth=1,  marker='o',alpha=0.9,edgecolors='black',
                    transform=ccrs.PlateCarree()
                    )
        else:
            bio.plot.scatter(x='longitude', y='latitude',s=4*bio['total_years'],ax=ax,label=b,color=colors2[n],
                    linewidth=1,  marker='o',alpha=0.9,edgecolors='black',
                    transform=ccrs.PlateCarree()
                    )
        # for i,j in bio.iterrows():
        #     x,y = plot_crs.transform_point(j.longitude, j.latitude, ccrs.PlateCarree())
        #     ax.annotate(str(j.biome)+str(j.total_years),(x,y))

    ax.set_extent([-180,180,-60,90],crs=ccrs.PlateCarree())

    plt.text(100,-50,str(str(numsites) + " sites"),fontsize=22,transform=ccrs.PlateCarree())

    #Improve the legend
    biomepatches=[]
    for bio in biomedict.keys():
        biomepatches.append(mpatches.Patch(color=colors2[biomedict[bio]], label=bio)) 
    legend = ax.legend(handles=biomepatches, title="Biomes", prop = {'size':'x-large'},
            handletextpad=1, columnspacing=1,
            loc="lower left", ncol=1, frameon=True, fontsize=22)

    labels = ["2", "5", "12", "20"]
    l1 = ax.scatter([],[], s=8, edgecolors='none',facecolors='k')
    l2 = ax.scatter([],[], s=20, edgecolors='none',facecolors='k')
    l3 = ax.scatter([],[], s=48, edgecolors='none',facecolors='k')
    l4 = ax.scatter([],[], s=80, edgecolors='none',facecolors='k')
    leg = ax.legend([l1, l2, l3, l4], labels, ncol=4, frameon=True, fontsize=16,
    handlelength=0.5, loc = 8, borderpad = 0.5, 
    handletextpad=0.3, scatterpoints = 1)  
    leg.set_title('Site Years',prop={'size':'x-large'})

    for item in ([ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(16)

    plt.setp(legend.get_title(),fontsize='18')

    # Add legends
    plt.gca().add_artist(legend)
    plt.gca().add_artist(leg)

    # # %%
    plt.savefig(str(outputFiles + 'mapofFLUXNETsites-world-biome.png'),dpi=450, bbox_inches='tight')
    plt.show()


def landscape_plot(df):
    plt.figure(figsize=(4*4, 4*5))

    plot_crs = ccrs.Robinson(central_longitude=0)
    ax = plt.axes(projection=plot_crs)
    ax.set_global()

    ax.coastlines('110m', linewidth=0.8, zorder=2,color='grey')

    for n,b in enumerate(landscapedict):
        lscape=df[df['landscape_state'].str.match(b)]
        lscape.plot.scatter(x='longitude', y='latitude',s=4*lscape['total_years'],ax=ax,label=b,color=colors3[n],
                linewidth=1,  marker='o',alpha=0.9,edgecolors='black',
                transform=ccrs.PlateCarree()
                )
    
    ax.set_extent([-180,180,-60,90],crs=ccrs.PlateCarree())

    # texts = [plt.text(df.longitude[i], df.latitude[i], df.Site[i], ha='center',
    #            va='center',transform=ccrs.PlateCarree()) for i in range(len(df.Site))]

    #ax2.set_boundary(circle, transform=ax2.transAxes)

    plt.text(100,-50,str(str(numsites) + " sites"),fontsize=22,transform=ccrs.PlateCarree())

    #Improve the legend
    lspatches=[]
    for lscape in landscapedict.keys():
        lspatches.append(mpatches.Patch(color=colors3[landscapedict[lscape]], label=lscape)) 
    legend = ax.legend(handles=lspatches, title="Landscape state", prop = {'size':'x-large'},
            handletextpad=1, columnspacing=1,
            loc="lower left", ncol=1, frameon=True, fontsize=22)

    labels = ["2", "5", "12", "20"]
    l1 = ax.scatter([],[], s=8, edgecolors='none',facecolors='k')
    l2 = ax.scatter([],[], s=20, edgecolors='none',facecolors='k')
    l3 = ax.scatter([],[], s=48, edgecolors='none',facecolors='k')
    l4 = ax.scatter([],[], s=80, edgecolors='none',facecolors='k')
    leg = ax.legend([l1, l2, l3, l4], labels, ncol=4, frameon=True, fontsize=16,
    handlelength=0.5, loc = 8, borderpad = 0.5, 
    handletextpad=0.3, scatterpoints = 1)  
    leg.set_title('Site Years',prop={'size':'x-large'})

    for item in ([ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(16)

    plt.setp(legend.get_title(),fontsize='18')

    # Add legends
    plt.gca().add_artist(legend)
    plt.gca().add_artist(leg)

    plt.savefig(str(outputFiles + 'mapofFLUXNETsites-world-landscapeStatus.png'),dpi=450, bbox_inches='tight')
    plt.show()


def peatland_plot(df):
    plt.figure(figsize=(4*4, 4*5))

    plot_crs = ccrs.Robinson(central_longitude=0)
    ax = plt.axes(projection=plot_crs)
    ax.set_global()

    ax.coastlines('110m', linewidth=0.8, zorder=2,color='grey')

    for n,b in enumerate(peatdict):
        ptlnd=df[df['peatland'].str.match(b)]
        ptlnd.plot.scatter(x='longitude', y='latitude',s=4*ptlnd['total_years'],ax=ax,label=b,color=colors3[n],
                linewidth=1,  marker='o',alpha=0.8,edgecolors='black',
                transform=ccrs.PlateCarree()
                )
    
    ax.set_extent([-180,180,-60,90],crs=ccrs.PlateCarree())

    # texts = [plt.text(df.longitude[i], df.latitude[i], df.Site[i], ha='center',
    #            va='center',transform=ccrs.PlateCarree()) for i in range(len(df.Site))]

    #ax2.set_boundary(circle, transform=ax2.transAxes)

    plt.text(100,-50,str(str(numsites) + " sites"),fontsize=22,transform=ccrs.PlateCarree())

    #Improve the legend
    ppatches=[]
    for n,b in enumerate(peatdict):
        ppatches.append(mpatches.Patch(color=colors3[n], label=b)) 
    legend = ax.legend(handles=ppatches, title="Peatland/Upland", prop = {'size':'x-large'},
            handletextpad=1, columnspacing=1,
            loc="lower left", ncol=1, frameon=True, fontsize=22)

    labels = ["2", "5", "12", "20"]
    l1 = ax.scatter([],[], s=8, edgecolors='none',facecolors='k')
    l2 = ax.scatter([],[], s=20, edgecolors='none',facecolors='k')
    l3 = ax.scatter([],[], s=48, edgecolors='none',facecolors='k')
    l4 = ax.scatter([],[], s=80, edgecolors='none',facecolors='k')
    leg = ax.legend([l1, l2, l3, l4], labels, ncol=4, frameon=True, fontsize=16,
    handlelength=0.5, loc = 8, borderpad = 0.5, 
    handletextpad=0.3, scatterpoints = 1)  
    leg.set_title('Site Years',prop={'size':'x-large'})

    for item in ([ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(16)

    plt.setp(legend.get_title(),fontsize='18')

    # Add legends
    plt.gca().add_artist(legend)
    plt.gca().add_artist(leg)

    plt.savefig(str(outputFiles + 'mapofFLUXNETsites-world-peatland.png'),dpi=450, bbox_inches='tight')
    plt.show()

### Generate all plots: ###
df = read_sites()

# Find the unique biomes we have and number of sites
biomes=df.biome.unique()
numsites= len(df)

biome_plot(df)
landscape_plot(df)
peatland_plot(df)
