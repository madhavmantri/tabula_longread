import matplotlib.pyplot as plt
import matplotlib.cm

# create a color dictionary for tissues
def tissue_colors():
    
    tissue_color_dict = {
             'Bladder': '#e7969c',
             'Blood': '#d6616b',
             'Bone_Marrow': '#cedb9c',
             'Ear': '#a1b84f',
             'Eye': '#c7ea46',#"#00ff7f",
             'Fat': '#e7cb94',
             'Heart': '#ff0800',
             'Kidney': '#7b4173',
             'Large_Intestine': '#31a354',
             'Liver': '#000080',
             'Lung': '#3182bd',
             'Lymph_Node': '#8c6d31',
             'Mammary':'#ce6dbd',
             'Muscle': '#e7ba52',
             'Ovary':'#fdc33c',
             'Pancreas': '#fd8d3c',
             'Prostate':'#637939',#'#a55194',#
             'Salivary_Gland':'#622a0f',
             'Skin': '#de9ed6',
             'Small_Intestine': '#6baed6',
             'Spleen': '#393b79',
             'Stomach':'#446478',
             'Testis':'#9cc4de',
             'Thymus': '#9c9ede',
             'Tongue':'#b5cf6b',
             'Trachea': '#969696',
             'Uterus':'#c64b8c',#'#ff0090',
             'Vasculature': '#843c39'
    }
    
    return tissue_color_dict




# create a color dictionary for donors
def donor_colors():
    donors = ['TSP1','TSP2','TSP3','TSP4','TSP5','TSP6','TSP7','TSP8','TSP9','TSP10','TSP11','TSP12','TSP13','TSP14','TSP15','TSP16','TSP17','TSP19','TSP20','TSP21','TSP25','TSP26','TSP27','TSP28','TSP30']
    
    import matplotlib.colors as pltcolors
    
    cmap = plt.get_cmap("YlGnBu")
        
    donor_color_dict = {}
    j=1/len(donors)
    for d in donors:
        donor_color_dict[d] = pltcolors.to_hex(cmap(j))
        j+=1/len(donors)
        
    return donor_color_dict


# create a color dictionary for donors
def compartment_colors():
    
    compartments = ['endothelial', 'epithelial', 'immune', 'stromal', 'germline','neural']
    
    import matplotlib.colors as pltcolors
    
    cmap = plt.get_cmap("YlOrRd")
        
    compartment_color_dict = {}
    j=1/len(compartments)
    for c in compartments:
        compartment_color_dict[c] = pltcolors.to_hex(cmap(j))
        j+=1/len(compartments)
        
    return compartment_color_dict


# create a color dictionary for methods
def method_colors():
    methods = ['10X','smartseq']
    
    import matplotlib.colors as pltcolors
    
    method_color_dict = {}
    method_color_dict['10X'] = '#90ee90'
    method_color_dict['smartseq2'] = '#006400'
    
    return method_color_dict

# create a color dictionary for methods
def assay_colors():
    assays = ['10X_5Prime_v2','10X_3Prime_v3.1','smartseq2','smartseq3']
    
    import matplotlib.colors as pltcolors
    
    assay_color_dict = {}
    assay_color_dict['10X_5Prime_v2'] = '#90ee90'
    assay_color_dict['10X_3Prime_v3.1'] = '#006400'
    assay_color_dict['smartseq2'] = '#b690ee'
    assay_color_dict['smartseq3'] = '#665480'
    
    return assay_color_dict

# create a color dictionary for sex
def sex_colors():
    sexes = ['male','female','undisclosed']
    
    import matplotlib.colors as pltcolors
    
    sex_color_dict = {}
    sex_color_dict['female'] = '#f4cae4'
    sex_color_dict['male'] = '#cbd5e8'
#     sex_color_dict['undisclosed'] = '#e6f5c9'
    
    return sex_color_dict

# create a color dictionary for SQANTI3 / pigeon structural categories
def structural_category_colors():

    # Canonical SQANTI3 category colors (FSM/ISM/NIC/NNC) plus the remaining
    # pigeon structural categories seen in this project's classification files.
    structural_category_color_dict = {
        'full-splice_match':       '#6BAED6',  # FSM - blue
        'incomplete-splice_match': '#FC8D59',  # ISM - orange
        'novel_in_catalog':        '#78C679',  # NIC - green
        'novel_not_in_catalog':    '#EE6A50',  # NNC - salmon/red
        'genic':                   '#969696',  # genic genomic - grey
        'genic_intron':            '#41B6C4',  # teal
        'antisense':               '#66C2A4',  # sea green
        'fusion':                  '#FFC125',  # goldenrod
        'intergenic':              '#E9967A',  # dark salmon
        'moreJunctions':           '#8DA0CB',  # muted blue-violet
    }

    return structural_category_color_dict


# visualize the color dictionaries
def plot_colortable(colors, title, sort_colors=True, emptycols=0):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    topmargin = 40

    # Sort colors by hue, saturation, value and name.
    by_hsv = [(v, k) for k, v in colors.items()]
    
    if sort_colors is True:
        by_hsv = sorted(by_hsv)
    names = [name for hsv, name in by_hsv]

    n = len(names)
    ncols = 4 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24, loc="left", pad=10)

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=colors[name], linewidth=18)

    return fig