# Classical Trajectories in Coulomb Three Body Systems

Authors:

- Fernando Perez <fernando.perez@berkeley.edu>
- Jorge Mahecha <mahecha@fisica.udea.edu.co>

Tools for the numerical study of classical Coulomb three body systems with
arbitrary masses and charges in three spatial dimensions.

The code contained in this repository can solve the electrostatic three body
problem (three charged particles) with a numerically stable scheme that uses a
regularized form of the equations of motion.

The code and manuscript contained in this repository were the senior
undergraduate thesis for F. Perez, conducted under the supervision of professor
Jorge Mahecha at the physics department of the Universidad de Antioquia, in
Medellin, Colombia.  This work was completed in 1994 and was part of the
following publications:

- F. Pérez and J. Mahecha. Classical trajectories in Coulomb three body
  systems. Rev. Mex. Física 42:6, 1070-1086 (1996). ([Full text PDF](
  http://rmf.smf.mx/pdf/rmf/42/6/42_6_1070.pdf)).

- A. Santander, J. Mahecha and F. Pérez, Rigid Rotator and Fixed shape
  solutions to the Coulomb Three-Body Problem. Few Body Systems, 22:1, 37-60
  (1997).
  
The `.tex` sources for the first manuscript, along with the original figures,
are included here for reference.

The purpose of this repository is to enable anyone interested in this research
to reproduce our results and further use this code.  We must note that this
repository contains unmodified files from 1994-1996, which are very unlikely to
compile today without additional work. But the binaries in the `c` directory
may work on a Windows machine or with a suitable DOS emulator.

Unfortunately, at this time the authors don't have any bandwidth to develop
this work any further.

The repository contains three directories:

- `c`: the main C sources for the programs, compiled at the time in Turbo
  C. Since these sources are unlikely to compile unmodified today, two binaries
  from that time are also included.  These can probably be made to run on a
  Windows computer with the right compatibility support (Windows XP most likely
  would work, I have no idea about newer versions).
  
- `maple`: the Maple worksheets for some of the analytical work, including the
  regularization of the equations of motion and the auto-generation of C form
  of the right-hand-side of these equations from the symbolic objects in Maple.
  
- `tex`: LaTeX sources for the manuscript submitted to the Revista Mexicana de
  Fisica, referenced above.


## License: BSD

The Coulomb Three Body code contained in this repository is released under the
terms of the Modified BSD License (also known as New or Revised or 3-Clause
BSD). See the accompanying COPYING file for details.
