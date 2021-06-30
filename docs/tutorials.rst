Tutorials
*********

Create customized oncoplots
===========================

We can use either the :command:`maf_oncoplt` command or the :meth:`pymaf.plot_oncoplot` method to create a "standard" oncoplot like the one shown below.

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/oncoplot.png

While the plot is pleasing to the eye, one may wish to customize it to add, for example, various annotation data for the samples like this:

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/customized_oncoplot_1.png

We can (relatively) easily achieve above thanks to the LEGO block-like plotting methods in the ``pymaf`` submodule (:download:`customized_oncoplot_1.py <examples/customized_oncoplot_1.py>`).

We can go one step further and sort the samples by one or more annotations (:download:`customized_oncoplot_2.py <examples/customized_oncoplot_2.py>`):

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/customized_oncoplot_2.png

Finally, we can also subset the samples with annotation data (:download:`customized_oncoplot_3.py <examples/customized_oncoplot_3.py>`):

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/customized_oncoplot_3.png

Control plot colors
===================

Let us consider the :meth:`fuc.api.pymaf.MafFrame.plot_snvclsc` method as an example. The method internally uses :meth:`seaborn.barplot` to draw bar plots.

.. plot::
    :context: close-figs

    from fuc import common, pymaf
    common.load_dataset('tcga-laml')
    maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
    mf = pymaf.MafFrame.from_file(maf_file)
    mf.plot_snvclsc()

In order to change the colors according to one of the `seaborn palettes <https://seaborn.pydata.org/generated/seaborn.color_palette.html#seaborn.color_palette>`__ (e.g. 'deep', 'muted', 'pastel'), we can use the ``palette`` option:

.. plot::
    :context: close-figs

    import seaborn as sns
    mf.plot_snvclsc(palette='pastel')

To choose one of the `colormaps available in maplotlib <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`__ (e.g. 'tab10', 'Dark2', 'Pastel1'), we can use the ``palette`` option with :meth:`seaborn.color_palette`:

.. plot::
    :context: close-figs

    mf.plot_snvclsc(palette=sns.color_palette('Pastel1'))

Some plotting methods use maplotlib instead of seaborn for coloring, like the :meth:`fuc.api.pymaf.MafFrame.plot_snvclss` method:

.. plot::
    :context: close-figs

    ax = mf.plot_snvclss(width=1)
    ax.legend(loc='upper right')

To choose a colormap from maplotlib:

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt
    ax = mf.plot_snvclss(width=1, color=plt.get_cmap('Pastel1').colors)
    ax.legend(loc='upper right')

To choose a palette from seaborn:

.. plot::
    :context: close-figs

    ax = mf.plot_snvclss(width=1, color=sns.color_palette('pastel'))
    ax.legend(loc='upper right')
