L. Zhang, K. Xiao, X. Wang, L. Kong, A novel fusion technology utilizing complex network and sequence information for FAD-binding site identification.*Analytical Biochemistry*, 2024, 285(15), 115401.  [https://doi.org/10.1016/j.ab.2023.115401](https://doi.org/10.1016/j.ab.2023.115401)

The benchmark datasets can be found in ./dataset, the codes for FADsite are available in ./src. And the PDB files of proteins are saved in ./pdb. Furthermore, the demo and corresponding documentation files can be found in ./demo. See our paper for more details.

Testing each proteins takes approximately 2 minute, depending on the sequence length.


### Testing the FADsite_seq on test4

```bash
cd ./src/
python FADsite_seq.py test4  
```
### Testing the FADsite_seq on test6

```bash
cd ./src/
python FADsite_seq.py test6 
```
### Test the FADsite on test4 (~6min)
```bash
cd ./src/
python FADsite.py test4  
```

### Testing the FADsite on test6 (~16min)
```bash
cd ./src/
python FADsite.py test6  
```
### contact
Kang Xiao: xiaokangneuq@163.com

