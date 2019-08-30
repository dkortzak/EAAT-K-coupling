# Allosteric gate modulation confers K<sup>+</sup>  coupling in glutamate transporters
Datasets and python scripts accompanying the paper *Allosteric gate modulation confers K<sup>+</sup>  coupling in glutamate transporters* by Kortzak et al.

### analysis scripts
  * [transition_path_eaat.py](../master/transition_path_eaat.py) reads discretized K<sup>+</sup> occupation state trajectories and calculates highest reactive flux pathways for EAAT1 (used in Fig. 2C). Depends on the [MSMexplorer](https://github.com/msmbuilder/msmexplorer) package. 
  * [transition_path_gltph.py](../master/transition_path_gltph.py) reads discretized K<sup>+</sup> occupation state trajectories and calculates highest reactive flux pathways for Glt<sub>Ph</sub> (used in Fig. 2D, E). Depends on the [MSMexplorer](https://github.com/msmbuilder/msmexplorer) package. 
  * [dwell_eaat1.py](../master/dwell_eaat1.py) reads discretized K<sup>+</sup> occupation state trajectories and calculates kinetic properties from EAAT1 simulations (used in Fig. 3A, B). 
  * [dwell_gltph_ofc.py](../master/dwell_gltph_ofc.py) reads discretized K<sup>+</sup> occupation state trajectories and calculates kinetic properties from Glt<sub>Ph</sub> OFC simulations (used in Fig. 3A, B).
  * [dwell_gltph_ifc.py](../master/dwell_gltph_ifc.py) reads discretized K<sup>+</sup> occupation state trajectories and calculates kinetic properties from Glt<sub>Ph</sub> IFC simulations (used in Fig. 3A, B).
  * [gen_states_ofc.py](../master/gen_states_ofc.py), [gen_states_ifc.py](../master/gen_states_ifc.py), [gen_states_eaat.py](../master/gen_states_eaat.py) describe the discretization procedure of distances into discretized K<sup>+</sup> occupation state trajectories.
  * [cgi_selectivity.py](../master/cgi_selectivity.py) analyzes data from alchemical transformations with the crooks-gaussian-intersection estimator (used in Fig. 3C). Script uses code copied from the [pmx package](https://github.com/dseeliger/pmx).
  * [cgi_mutant.py](../master/cgi_mutant.py) analyzes data from alchemical transformations with the crooks-gaussian-intersection estimator (used in Fig. 4A). Script uses code copied from [pmx package](https://github.com/dseeliger/pmx).
  * [eaat_profiles.py](../master/eaat_profiles.py) reads and plots probability profiles. Calculates and plots HP2 closed probablility for EAAT1 (used in Fig. 6 and 7).
  * [gltph_profiles.py](../master/gltph_profiles.py) reads and plots probability profiles. Calculates and plots HP2 closed probablility for Glt<sub>Ph</sub> (used in Appendix Fig. S7C).

### structure
  * PDB files to generate the images in Fig. 5C and Fig. 7A.
      * [EAAT1_apostruc.pdb](../master/structure/EAAT1_apostruc.pdb), [EAAT1_K1struc.pdb](../master/structure/EAAT1_K1struc.pdb) representative apo and K1 bound EAAT1 structures (used in Fig. 5).
      * [EAAT1_hpclosed_nosb.pdb](../master/structure/EAAT1_hpclosed_nosb.pdb), [EAAT1_hpopen_sb.pdb](../master/structure/EAAT1_hpopen_sb.pdb) representative EAAT1 structures of the HP2 open with saltbridge and HP2 closed without saltbridge (used in Fig. 7).
      
### data
* [freeMD_EAAT1](../master/data/freeMD_EAAT1), [freeMD_GLTPH_OFC](../master/data/freeMD_GLTPH_OFC), [freeMD_GLTPH_IFC](../master/data/freeMD_GLTPH_IFC) 
contain discretized K<sup>+</sup> occupation state trajectories from independent MD simulations (used in Fig. 2C, D, E).
* [cgi/selectivity](../master/data/cgi/selectivity)
contains work values from alchemical transformations of a bound ion for each binding site in Glt<sub>Ph</sub> and EAAT1 (used in Fig. 3C).
* [cgi/mutants](../master/data/cgi/mutants) contains work values from alchemical transformations of side-chains for apo and K1-bound Glt<sub>Ph</sub> (used in Fig. 4A).
* [umbrella_sampling_GLTPH](../master/data/umbrella_sampling_GLTPH) contains probability profiles for HP2 opening for OFC Glt<sub>Ph</sub>. 
* [umbrella_sampling_EAAT1](../master/data/umbrella_sampling_EAAT1) contains probability profiles for HP2 opening for EAAT1 WT, R479A, L448A, L448T EAAT1 in apo and K1-bound states (used in Fig. 6 and Fig. 7).

### citation
If you make use of this data, please cite it.
```
@article{kortzak2019,
  title={Principles of K+ coupling in the glutamate transporter family},
  author={Daniel Kortzak, Claudia Alleva, Ingo Weyand, David Ewers, Meike I. Zimmermann, Arne Franzen, Jan-Philipp Machtens and Christoph Fahlke},
  journal={...},
  pages={...},
  year={...},
  publisher={...}
}
