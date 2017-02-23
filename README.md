# WARNING! 
k-point restart results are __NOT CONSISTENT__ with CPMD and abinit!

Modified source code of QuantumESPRESSO 5.3.0 for alchemy purpose.
The detailed comparison can be visualized by github 
[compare branch](https://github.com/SamKChang/espresso-5.3.0-alchemy/compare/original-espresso-5.3.0...master).

Reference calculation: 
To save restart files for reference run, use 'ALCHEMY' directive with keyword 'reference', 
i.e. add "ALCHEMY reference" to pw.x input file

Alchemy prediction: 
To restart calculations for alchemy prediction, use 'ALCHEMY' directive with keyword 'prediction', 
i.e. add "ALCHEMY prediction" to pw.x input file
This requires restart files from reference run. 
_Note_: full k-point mesh is necessary if the symmetry of reference calculation and prediction calculation
are different.

Original README
===============

This is the distribution of the Quantum ESPRESSO suite of codes (ESPRESSO: 
opEn-Source Package for Research in Electronic Structure, Simulation, 
and Optimization), promoted by the IOM-DEMOCRITOS National Simulation Center 
of the Italian CNR (http://www.democritos.it). 

Quick installation instructions for the impatient:
   ./configure [options]
   make all
("make" alone prints a list of acceptable targets). Binaries go in bin/.
For more information, see the general documentation in directory Doc/, 
package-specific documentation in */Doc/, and the web site
http://www.quantum-espresso.org/

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.


