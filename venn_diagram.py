from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles


# Subset sizes
s = (
    50,    # Abc # UK
    200,    # aBc #DK
    300,    # ABc #UK/FR
    200,    # abC #FR
    400,    # AbC #DK/FR
    500,  # aBC #DK/FR/UK
    50,    # ABC #UK/FR
)

v = venn3(subsets=s, set_labels=('A', 'B', 'C'))

# Subset labels
v.get_label_by_id('100').set_text('UK')
v.get_label_by_id('010').set_text('DK')
v.get_label_by_id('110').set_text('UK-DK')
v.get_label_by_id('001').set_text('FR')
v.get_label_by_id('101').set_text('UK-FR')
v.get_label_by_id('011').set_text('DK-FR')
v.get_label_by_id('111').set_text('DK-FR-UK')

# Subset colors
v.get_patch_by_id('100').set_color('c')
v.get_patch_by_id('010').set_color('#993333')
v.get_patch_by_id('110').set_color('blue')

# Subset alphas
v.get_patch_by_id('101').set_alpha(0.4)
v.get_patch_by_id('011').set_alpha(1.0)
v.get_patch_by_id('111').set_alpha(0.7)

# Border styles
c = venn3_circles(subsets=s, linestyle='solid')
c[0].set_ls('dotted')  # Line style
c[1].set_ls('dashed')
c[2].set_lw(1.0)       # Line width

plt.show()