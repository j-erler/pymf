PyMF
=========

**PyMF** is a pure-Python implementation of several spatial filtering 
image processing techniques, namely the 
`matched filter <http://dx.doi.org/10.1093/mnras/279.2.545>`_, 
`matched multifilter <http://dx.doi.org/10.1046/j.1365-8711.2002.05704.x>`_ and the newly 
introduced `constrained matched filter <https://arxiv.org/abs/1809.06446>`_ and
`constrained matched mutifilter <https://arxiv.org/abs/1809.06446>`_. These techniques
are designed for applications in astronomy and use the flat-sky approximation.

This documentation provides an overview on the installation process and the available 
functions that are part of PyMF. Practical examples are provided in two jupyter 
notebooks that can be found in the /examples directory.

PyMF is being actively developed on `GitHub <https://github.com/j-erler/pymf>`_.

.. image:: https://img.shields.io/badge/GitHub-j--erler%2Fpymf-blue.svg?style=flat
    :target: https://github.com/j-erler/pymf
.. image:: https://img.shields.io/badge/docs-passing-green.svg?style=flat
    :target: https://pymf.readthedocs.io/en/latest/index.html#
.. image:: https://img.shields.io/badge/license-MIT-red.svg?style=flat
    :target: https://github.com/j-erler/pymf/blob/master/LICENSE    
.. image:: https://img.shields.io/badge/arXiv%3A-1809.06446-orange.svg?style=flat
    :target: https://arxiv.org/abs/1809.06446
    
.. toctree::
   :maxdepth: 2
   :caption: Contents:

Acknowledgement
---------------

Please cite `Erler Ramos-Ceja, Basu & Bertoldi (2018)
<https://arxiv.org/abs/1809.06446>`_ if you find this code useful in your
research.
The BibTeX entry for the paper is::

    @ARTICLE{2018arXiv180906446E,
        author = {{Erler}, J. and {Ramos-Ceja}, M.~E. and {Basu}, K. and {Bertoldi}, F.},
         title = "{Introducing constrained matched filters for improved separation of point sources from galaxy clusters}",
       journal = {ArXiv e-prints},
          year = 2018,
        eprint = {1809.06446}
    }

Installation
---------

.. toctree::
   :maxdepth: 2

   install


Reference
---------

.. toctree::
   :maxdepth: 2

   code

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

License
-------

.. toctree::
   :maxdepth: 1

   license
