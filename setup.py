from setuptools import setup, find_packages

setup(
  name = 'glycanPRMQuant',
  packages = find_packages(),
  version = '0.1.0',
  license='MIT',
  install_requires=[
    'numpy',
    'pandas',
    'scipy',
    'matplotlib',
    'seaborn',
    'statsmodels',
    'scikit-learn',
    'openpyxl',
    'scienceplots',
    'pyteomics',
    'glypy'
  ],
  entry_points={
    'console_scripts': [
      'glycan-prmquant=glycanPRMQuant.cli:main',
    ],
  },
  description = 'A package for glycan PRM quantification',
  author = 'Vishal Sandilya',
  author_email = 'vishal.sandilya@ttu.edu',
  python_requires='>=3.12'
)
