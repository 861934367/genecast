genecast
======

A command-line toolkit and Python library for scientific research cooperation

Support
=======

Please use email to ask any questions to zhou.tao@genecast.com.cn

Installation
============

genecast runs on >= Python 3.5. Your operating system might already provide Python
3.5, which you can check on the command line::

python --version

python setup.py install


Usage

============

genecast.py -h


Python dependencies
-------------------

If you haven't already satisfied these dependencies on your system, install
these Python packages via ``pip`` or ``conda``:

- `Biopython <http://biopython.org/wiki/Main_Page>`_
- `Reportlab <https://bitbucket.org/rptlab/reportlab>`_
- `matplotlib <http://matplotlib.org>`_
- `NumPy <http://www.numpy.org/>`_
- `SciPy <http://www.scipy.org/>`_
- `Pandas <http://pandas.pydata.org/>`_
- `sklearn`

On Ubuntu or Debian Linux::

    sudo apt-get install python-numpy python-scipy python-matplotlib python-reportlab python-pandas
    sudo pip install biopython pyfaidx sklearn pyvcf --upgrade

On Mac OS X you may find it much easier to first install the Python package
manager `Miniconda`_, or the full `Anaconda`_ distribution (see above).
Then install the rest of CNVkit's dependencies::

    conda install numpy scipy pandas matplotlib reportlab biopython pyfaidx sklearn

    pip install numpy scipy pandas matplotlib reportlab biopython pyfaidx sklearn




