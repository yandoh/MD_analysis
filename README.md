# MD_analysis
analysis codes for MD trajectory generated by modylas program.
(about modylas, see https://www.modylas.org/)

1. atomgr      ... A Fortran code to calculate g(r) between atom pairs.
2. copys       ... A Fortran code to copy centers of mass of the all molecules in a calculation unit cell into its image cells.
3. dcdmol      ... A Fortran code to recover shape of each molecule in .dcd file.
4. mdmntrave   ... A Fortran code to take an averave of values in .mdmntr files.
5. molgr       ... A Fortran code to calculate g(r) between molecules' mass centers.
6. partgr       ... A Fortran code to calculate g(r) between mass centers for parts (segments) of molecules. [A part can be defined arbitrarily by a user]
7. SelfD       ... A Fortran code to calculate the self-diffusion coefficient of molecules.
8. OrderPs     ... A Fortran code to calculate the order parameter to evaluate stacking of two plane segments (e.g., benzene rings)
9. viscosity   ... A Fortran code to calculate the viscosity by two methods (viscosity1 and viscosity2) based on pressure tensor.
10. catdcd      ... A link for download cite of the catdcd program.
