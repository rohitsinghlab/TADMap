from setuptools import setup

setup(name='tadmap',
      version='0.1.2',
      description='Consensus layout of topologically associating domains (TADs) and utilities for RNA-seq analysis',
      
      long_description='Topologically associating domains (TADs) are contiguous segments of the genome where the genomic elements are in frequent contact with each other. The TAD Map provides a consensus estimate, aggregated from multiple experimental datasets, of this layout in human and mouse. It also provides tools to map any single-cell RNA-seq dataset to TAD signatures, with gene expression mapped to TAD activation probabilities in each cell',      
      url='http://github.com/rs239/tadmap',
      author='Rohit Singh',
      author_email='rsingh@alum.mit.edu',
      license='MIT',
      packages=['tadmap', 'tadmap.data'],
      install_requires = 'numpy,scipy,pandas,sklearn,tqdm,scanpy'.split(','),
      zip_safe=False)
