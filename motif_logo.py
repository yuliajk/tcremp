
import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
fp = FontProperties(family="Arial", weight="bold") 
globscale = 1.35

LETTERS = { "A" : TextPath((-0.305, 0), "A", size=1, prop=fp),
            "R" : TextPath((-0.384, 0), "R", size=1, prop=fp),
           "N" : TextPath((-0.384, 0), "N", size=1, prop=fp),
           "D" : TextPath((-0.384, 0), "D", size=1, prop=fp),
           "C" : TextPath((-0.384, 0), "C", size=1, prop=fp),
           "Q" : TextPath((-0.384, 0), "Q", size=1, prop=fp),
           "E" : TextPath((-0.384, 0), "E", size=1, prop=fp),
           "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
           "H" : TextPath((-0.384, 0), "H", size=1, prop=fp),
           "I" : TextPath((-0.384, 0), "I", size=1, prop=fp),
           "L" : TextPath((-0.384, 0), "L", size=1, prop=fp),
           "K" : TextPath((-0.384, 0), "K", size=1, prop=fp),
            "M" : TextPath((-0.35, 0), "M", size=1, prop=fp),
           "F" : TextPath((-0.35, 0), "F", size=1, prop=fp),
           "P" : TextPath((-0.35, 0), "P", size=1, prop=fp),
           "S" : TextPath((-0.35, 0), "S", size=1, prop=fp),
           "T" : TextPath((-0.35, 0), "T", size=1, prop=fp),
           "W" : TextPath((-0.35, 0), "W", size=1, prop=fp),
           "Y" : TextPath((-0.35, 0), "Y", size=1, prop=fp),
           "V" : TextPath((-0.35, 0), "V", size=1, prop=fp),
           "-" : TextPath((-0.35, 0), "-", size=1, prop=fp)}

LETTERS_color =  {"A" : 'black',
            "R" : 'blue',
           "N" : 'purple',
           "D" : 'orange',
           "C" : 'green',
           "Q" : 'purple',
           "E" : 'orange',
           "G" : 'green',
           "H" : 'blue',
           "I" : 'black',
           "L" : 'black',
           "K" : 'blue',
            "M" : 'black',
           "F" : 'black',
           "P" : 'black',
           "S" : 'green',
           "T" : 'green',
           "W" : 'black',
           "Y" : 'green',
           "V" : 'black',
           "-" : 'black'}

def amino_letterAt(letter, x, y, yscale=1, ax=None):

    text = LETTERS[letter]

    t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
        mpl.transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=LETTERS_color[letter],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p


def plot_amino_logo(motif, title, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,3))
    
    x = 1
    maxi = 0
    
    for xi in range(motif.shape[1]):
        scores = motif.iloc[:, xi]
        y = 0
        
        for base,score in scores.items():
            if score > 2/len(scores):
                amino_letterAt(base, x,y, score, ax)
                y += score
        x += 1
        maxi = max(maxi, y)
    #plt.xticks(range(1,x))
    ax.set_xlim((0, x)) 
    #plt.ylim((0, maxi)) 
    #plt.suptitle(title)
    #plt.tight_layout()      
    #plt.show()