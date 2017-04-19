
from distutils.core import setup

from os.path import dirname, join
from glob import glob
import sh

setup_args = {}
DIR = (dirname(__file__) or '.')
from distutils.core import setup
setup_args.update(
    name='genecast',
    version='0.1.2',
    author='zhou tao',
    author_email='zhou.tao@genecast.com.cn',
    packages=[
    'genecast_package',
    "genecast_package.sklearn"],
    scripts=[join(DIR, 'genecast.py')] + glob(join(DIR, 'scripts/*.py')))

setup(**setup_args)
