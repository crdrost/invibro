# invibro

This program was used to generate several pretty graphs for my Master's thesis,
"Vibration-Dependent Quantum Transport". It is based on an approach previously
derived by Karsten Flensberg in "Tunneling broadening of vibrational sidebands
in molecular transistors," *Phys. Rev. B* 68: 205323.

The file structure is a simple Python module; it is intended to run on Python
2.7.3 with access to the Python packages for `scipy` (and therefore `numpy`) and 
`Gnuplot` -- the latter is used only by the examples, and can be replaced with
other packages (I made significant use of `matplotlib` when browsing around 
phase space at low resolutions). 

This means that under Ubuntu the dependencies can probably be entirely installed
via:

    sudo apt-get install python python-gnuplot gnuplot python-scipy

The examples assume a Unix environment, namely a `/tmp` folder into which 
results are printed as postscript files.

# License

This project was authored by Chris Drost of drostie.org. To the extent 
possible by all laws in all countries, I hereby waive all copyright and any 
related rights under the Creative Commons Zero (CC0) waiver/license, which 
you may read online at:

    http://creativecommons.org/publicdomain/zero/1.0/legalcode

This means that you may copy, distribute, modify, and use my code without 
any fear of lawsuits from me. It also means that my code is provided with NO
WARRANTIES of any kind, so that I may have no fear of lawsuits from you. 

This waiver/license applies to all of the code in the project as well as this 
particular file. 
