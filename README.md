# Not supported

Please not I do not maintain this code any longer, nor do I have access to the `optics_5.4.1` code. If you do not already have that code then it is unlikely you can use this tool.

I recommend the `sl3me` tool which does the same but is supported: https://github.com/ldwillia/SL3ME

# SLMETools
The SLME method is outlined in the paper of Yu and Zunger, [PRL, 108,068701](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.108.068701)

The following is a collection of some of the tools you need to calculate the SLME.

* Run a normal VASP calculation using INCAR.ele and dense, $$\Gamma$$-point centred k-mesh.
* Obtain the minimum band gap [we call this E_g^fund]
* Copy IBZKPTS to KPOINTS
* Run a VASP calculation using INCAR.opt [Produces OPTIC]
* Copy the OPTCTR.eps file to OPTCTR
* Run the optics_5.4.1 [Produces EPS]
* Run the script make_diel.py [Produces DIEL]
* Copy the OPTCTR.trans to OPTCTR
* Run `optics_5.4.1 > Transitions.dat` [Returns the value of E_g^allowed]
* Now run the SLME code: `python calc_slme.py DIEL E_g^fund E_g^allowed 0.1`


NB the optics processing is available from the [VASP site](http://www.freeware.vasp.de/VASP/optics/)
