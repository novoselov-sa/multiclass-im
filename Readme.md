This is a code for class group computation for imaginary multiquadradic fields.

It is built on top of the implementation of the Biasse-van Vredendaal algorithm for the class group computation for real multiquadratic fields:
* https://scarecryptow.org/publications/multiclass.html

# Requirements
* Sage 9.4+

# How to compile
1. Run ```make```
2. Edit paths in Makefile according to your installation if code above fails.

# Folder structure
* tests/ - unit tests
* trees/ - trees describing prime ideal splitting over subfields. Factor base is loaded from these files. They are generated by ```trees_generation.sage``` script.
* relations/ - class group relation matrices for every subfield. They are generated by ```testrelations.sage``` script.
* experiments/ - examples of class group computation.

# How to use
1. Set the target field in the file ```trees_generation.sage``` (```food``` variable). Alternatively, you can pass d_1 d_2 ... d_n as command line arguments of this script.
2. Run ```sage trees_generation.sage``` to generate trees.
3. Run ```sage testrelations.sage``` to compute class group and relation matrices for subfields.
4. Run ```sage clgp_verify.sage d_1 d_2 ... d_n --log-zeta-res z --class-number h``` to verify result of computation. Here z is the residue of Dedekind zeta function at 1. This value can be computed using [Hecke](https://www.thofma.com/Hecke.jl/dev/), see examples in ''experiments'' folder.