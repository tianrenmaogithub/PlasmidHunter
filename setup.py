import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        try:
            subprocess.check_call(['conda', 'install', '-c', 'bioconda', 'diamond=2.1.8'])
            subprocess.check_call(['conda', 'install', '-c', 'bioconda', 'prodigal'])
        except subprocess.CalledProcessError:
            print("Failed to install the required dependencies (diamond ==2.1.8, prodigal) using Conda. Please install them manually.")

setup(
    name='plasmidhunter',
    version='1.3',
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
    package_data={'': ['model/*pkl']},
    entry_points={
        'console_scripts': [
            'plasmidhunter=PlasmidHunter.PlasmidHunter:main',
        ],
    },
    cmdclass={
        'install': PostInstallCommand,
    }

)

