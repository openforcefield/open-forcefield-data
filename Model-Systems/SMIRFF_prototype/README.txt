Please see smarff_mockup1.Frosst_AlkEthOH.txt for more extensive documentation.

This README is from an e-mail from Christopher Bayly describing contents of this directory:

Please find attached the first SMIRFF ever (actually a SMARFF 'cuz I only use smarts not smirks yet) and some associated useful files. It is only a mockup because I have no infrastructure yet with which to even parse the file to see if it is well-formed.

I think this is directed particularly to John so we can dialog over the format; I think John will have to change it and that is fine; this is trying to get the ball rolling to get the SMIRFF infrastructure going.

The SMARFF is smarff_mockup1.Frosst_AlkEthOH with important documentation in smarff_mockup1.Frosst_AlkEthOH.txt.

To help understand the smarts strings I have also attached the .py dump of my iPython notebook smartsDepict.ipynb (also attached) which takes a given smarts and highlights wherever it appears in the structure... it uses OEDepict in the stream so that I can see it in iPython but I don't think it will work as a standalone script. Its what I have so I am giving it to you to work from. I used it with my "contains all parameters at least once" test molecule AlkEthOH_dvrs1.oeb (not in the toy dataset but maybe should be) which I am also attaching.

Beyond this I want to do two other mockups but I doubt I can get to it this week:
- smarff_mockup2.simple_AlkEthOH: a greatly simplified version for AlkEthOH set, possibly as a starting point for automated smarts modification
- smirff_mockup1.simple_AlkEthOH: a SMIRFF version of smarff_mockup2.simple_AlkEthOH, beginning to use the atom-mapping feature of SMIRKS, possibly as a starting point for automated smarts modification

