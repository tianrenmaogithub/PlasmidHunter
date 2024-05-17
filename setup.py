import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='plasmidhunter',
    version='1.4',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=[
        'requests',
        'pandas',
        'numpy==1.25.0',
        'biopython',
        'scikit-learn==1.3.2',
    ],
    python_requires='==3.10.*',    
    include_package_data=True,
    package_data={'': ['model/*gz']},
    entry_points={
        'console_scripts': [
            'plasmidhunter=PlasmidHunter.PlasmidHunter:main',
        ],
    }
)

