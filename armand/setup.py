# Setup.py python distribution script to install the NOM DEL PROJECTE in the computer

from distutils.core import setup

setup(
   name='NOM DEL PROJECTE',
   version='1.0',
   description='NOM_DEL_PROJECTE is a python program designed to generate macrocomplex structures from simple pair interactions',
   long_description=open('README.md').read(),
   author='Oscar Camacho, Armand González',
   author_email = 'osc1997@gmail.com, armand.gonzalez01@estudiant.upf.edu',
   url='https://github.com/',
   packages=['MB'],
   install_requires=['biopython >= 1.73.0',
		                 'numpy >= 1.16.1'],
   license='LICENSE.txt',
   classifiers=[
                "Programming Language :: Python :: 3",
	              "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent"],
    scripts=["MB/MBlauncher.py" , "MB/MB_GUI.py"]
)


from distutils.core import setup

setup(name='BYCProject', version='1.0', description='Build a biological macrocomplex from pairwise interactions',
      author='Miquel Angel Schikora Tamarit and Marina Lleal Custey',
      author_email='miquelangel.schikora01@estudiant.upf.edu, marina.lleal01@estudiant.upf.edu',
      url='https://github.com/MikiSchikora/SBI_project',
      packages=['byc'], scripts=['byc/run_byc.py','byc/refine_model.py'], requires=['Bio'])
