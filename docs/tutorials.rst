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
