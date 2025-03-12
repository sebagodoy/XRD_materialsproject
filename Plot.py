import os
import json
import matplotlib.pyplot as plt
import numpy as np

# Peaks as Gauss-Lorentz mixing (only esthetic, not formal Rietveld)
MixFact = .5        # mixing factor f*Gauss + (1-f)*Lorentz
WidthGauss = 1      # width Gauss peak
WidthLorentz = .8   # width Lorentz peak

# Quality of the plotting
extend = 5          # extend 2theta range of ploting
spacing = 10000     # resolution

# Include annotations
ShowPeakNumber = True   # Include peak position
ShowPlanes = True       # Include plane definition

# ------------------------------------------------------------------------------------------------------------------

def Lorentz(_th, _thcenter, Amp=100):
    gamma = WidthLorentz
    return (Amp*np.pi*gamma)*(1/(gamma * np.pi)) /(1+((_th - _thcenter)/gamma)**2)

def Gauss(_th, _thcenter, Amp=100):
    cfact = WidthGauss
    return Amp * np.exp(-((_th - _thcenter)**2)/(2*cfact**2))

def GaussLorentz(_th, _thcenter, Amp=100, Mix=.5):
    return Mix * Gauss(_th, _thcenter, Amp) + (1-Mix) * Lorentz(_th, _thcenter, Amp)

FileList = [i for i in os.listdir('./') if '.json' in i]
FileList.sort()
print(FileList)

fig,axs = plt.subplots(len(FileList),1, sharex=True)
plt.subplots_adjust(hspace=0, top=.98, right=.98, left=.02)


for iF, iFile in enumerate(FileList):
    plt.sca(axs[iF])
    axs[iF].set_ylim(0, 100+extend)

    # Open and parse
    with open(iFile, 'r') as file:
        data = json.load(file)
    _theta = [i[2] for i in data['pattern']]
    _Ampl = [i[0] for i in data['pattern']]

    # Title
    if "!" in iFile.split('json')[0]:
        plt.annotate(iFile.split('json')[0].split('!')[0], xy=(.5, .98), xycoords='axes fraction',
                     ha='center', va='top',
                     bbox=dict(facecolor='white', edgecolor='none', alpha=.5))

        plt.annotate(iFile.split('json')[0].split('!')[1],
                     xy=(.98, .98), xycoords='axes fraction',
                     ha='right', va='top', fontsize='x-small',
                     bbox=dict(facecolor='white', edgecolor='none', alpha=.5))
    else:
        plt.annotate(iFile.split('json')[0], xy=(.5, .98), xycoords='axes fraction',
                     ha='center', va='top',
                     bbox=dict(facecolor='white', edgecolor='none', alpha=.5))

    plt.annotate(f"({data['wavelength']['element']} , {data['wavelength']['in_angstroms']}"
                 + r"$\AA$)",
                 xy=(.02, .98), xycoords='axes fraction',
                 ha='left', va='top', fontsize='x-small',
                 bbox=dict(facecolor='white', edgecolor='none', alpha=.5))

    # Patterns
    for i in data['pattern']:
        plt.axvline(x=i[2], ymin=0, ymax=i[0]/(100+extend), linewidth=.8, linestyle='solid')
        plt.axvline(x=i[2], ymin=0, c='k', alpha=.5, linewidth=.5, linestyle='dashed')

        if ShowPeakNumber or ShowPlanes:
            txt = ' '
            if ShowPeakNumber:
                txt += '{:.1f}'.format(i[2])
            if ShowPeakNumber and ShowPlanes:
                txt += ' , '
            if ShowPlanes:
                txt += str(i[1]) + ', d='+'{:.2f}'.format(i[3])
            # identification
            plt.annotate(txt,
                         xy=(i[2], .0),
                         xycoords=('data', 'axes fraction'),
                         rotation=90, fontsize=5, ha='right', va='bottom')

    # Shade
    _theta = [i[2] for i in data['pattern']]
    _Ampl = [i[0] for i in data['pattern']]
    xlimsXRD = (min(_theta)-5, max(_theta)+5)
    _xXRD = [xlimsXRD[0] + i*(xlimsXRD[1]-xlimsXRD[0])/spacing for i in range(spacing)]


    _yXRD = [sum([GaussLorentz(i, j, k, MixFact) for j,k in zip(_theta, _Ampl)]) for i in _xXRD]


    plt.fill_between(_xXRD, _yXRD, color='r', alpha=0.2)

    axs[iF].get_yaxis().set_visible(False)



plt.gca().invert_xaxis()
plt.xlabel(r'$2\theta$')


plt.show()