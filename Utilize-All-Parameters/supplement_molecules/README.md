## Adding supplemental molecules

In this directory, I (David Mobley) an adding supplemental molecules to ensure we cover all parameters, based on finding specific molecules which utilize parameters which were not covered in the first pass.

I am adding supplemental molecules noted below to the `chosen.smi` from ../selected, as of June 27, 2019 after merged of the updated `chosen.smi`..

### Coverage testing:
- To test coverage, run `python ../scripts/check_param_coverage.py -n N -t N --start 0 --end N-1 -f chosen_supplemented.smi -d .` where N is the number of molecules in the `chosen_supplemented.smi` file; this outputs the parameter usage for the set.

### Molecules added/issues fixed:
 - `b47`, for a sulfur-iodine bond, was not utilized so we added `Cc1ccc(cc1)S(=O)(=O)I`, an enamine building block/eMolecules 50366321, 4-methylbenzene-1-sulfonyl iodide.
- `b75`, for nitrogen-fluorine bonds, was not utilized so we added `c1ccc(c(c1)S(=O)(=O)N)S(=O)(=O)NF`, N-fluoro-o-benzenedisulfonamide. 
- `b78`, a nitrogen-iodine bond, was not utilized so we added `C1CC(=O)N(C1=O)I`, N-iodosuccinimide.
- `b82`, a phorphorus-iodine bond, was not utilized so we added `c1ccc(cc1)P(c2ccccc2)(c3ccccc3)(I)I`, eMolecules 506528, Diiodotriphenylphosphorane.
- `t82`, an azide torsion, was not utilized so we added `c1ccc(cc1)CN=[N+]=[N-]`, azidomethylbenzene .
- `t182a`, a tetrazole torsion added, was not utilized so added the tetrazole in question, `Cc1ccccc1c2[n-]nnn2`, 5-(o-tolyl)-1H-tetrazole.
- `b37`, N#N bond, added C1=NC(C(=N1)C(=O)N)[N+]#N (emolecules 45917281)
- `b51`: `[#16X2:1]-[#8X2:2]`. Somewhat unusual but adding methylsulfanol, CSO.
- b61: `[#15:1]=[#7:2]`. Not very common. Adding c1ccc(cc1)P(=N)(c2ccccc2)c3ccccc3, eMolecules 32261273
- b81: `[#15:1]-[#35:2]`. Not very common. Add tribromophosphane (P(Br)(Br)Br)
- t155: `[#6X3:1]-[#7:2]-[#15:3]=[*:4]`. Add orthene, CC(=O)NP(=O)(OC)SC, eMolecules 502047

That's 82 lines, and these additions do remove those instances as missing parameters. 
