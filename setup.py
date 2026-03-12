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
    'networkx'
  ],
  extras_require={
    'gui': ['PyQt5>=5.15.0'],
  },
  entry_points={
    'console_scripts': [
      'glycan-builder-gui=glycanPRMQuant.glycanBuilderGUI:launch_gui',
    ],
  },
  description = 'A package for glycan PRM quantification',
  author = 'Vishal Sandilya',
  author_email = 'vishal.sandilya@ttu.edu',
  python_requires='>=3.12'
)