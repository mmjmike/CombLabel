# UUCSL Package

Implemented in Python 3, **scipy** package required

Install Python:
https://www.python.org/downloads/

How to install python packages:
https://packaging.python.org/tutorials/installing-packages/


## Contributors

1. Mikhail Yu. Myshkin <mikhail.myshkin@phystech.edu>
2. Maxim A. Dubinnyi <dumkas@yandex.ru>
3. Zakhar O. Shenkarev <zakhar-shenkarev@yandex.ru>

Group of Structural Biology of the Ion Channels, IBCh RAS, Russia,
Moscow <http://www.ibch.ru/en/structure/groups/sbic>

---
## Licence

The program is distributed under AGPL-3.0 license and is
available at https://github.com/mmjmike/UUCSL


---
## Goal of the program

UUCSL is a Python 3 package that calculates universal unambiguous
combinatorial selective labeling schemes to label proteins with stable
13C and 15N isotopes and perform fast NMR signal assignment

---

## Package scripts

* **ncsgen.py** - interactive console application, that generates files,
    describing NMR coding systems, based on the spectra and labeling types
    which are entered by user
* **blockfinder.py** - generates Elementary Blocks (ELB) for universal
    scheme calculation, based on the chosen NCS and number of samples limit
* **elbclean.py** - cleans redundant blocks from elementary blocks file
* **elbstats.py** - calculates stats for elementary blocks files
* **UUCSL.py** - main script that calculates universal labeling schemes
    for a given list of amino acid types and ELB file, and performs
    price optimization of the schemes based on amino acid prices and
    availability of particular labeling types
* **CombLabel_check.py** - calculates labeling dictionary for particular
    protein sequence and UUCSL solution and prints out solution statistics

---

## How to use

### ncsgen.py

**ncsgen.py** is run without any command line arguments and will
inetractively ask
the user to enter particular data about NMR coding system. After
obtaining all required information it automatically generates `*.ncs`
file with following contents:

~~~~
[NCS = NC2]
[Deuterated = False]
[Labels = X,N,C]
[Spectra = HSQC,HNCO]

[code_pairs]
 , X, N, C
X, 0, 1, 0
N, 0, 1, 0
C, 0, 2, 0

[codes]
Code,   HSQC,   HNCO
   0,      0,      0
   1,      1,      0
   2,      1,      1
~~~~
Typical NCSs are collected in `NCS` folder of the package.


---
### blockfinder.py

**blocksfinder.py** searches for ELementary Blocks (ELB), out of which
 the labeling scheme can be constructed. As the input it takes two
 positional command line arguments: `ncs` (NCS) and `samples` (maximal
 number of
 samples in the resulting blocks). It starts from one-sample blocks and
 moves to more samples as it finds all blocks in previous number of
 samples and outputs non-equivalent blocks with at least two patterns
 to file `*_elb.txt`. It is wise to run this script with low number of
 samples first (2-4) before increasing the number, as the task of
 finding the blocks becomes combinatorially hard.

 **blocksfinder.py** has some more command-line arguments:

 * `--minpatterns, -p` - minimal number of patterns that outputted ELBs
  must have. Slightly increases the speed of the calculation and is
  suitable for checking whether a labeling scheme of P patterns can be
  found in S samples for a given NCS
 * `--verbose, -v` - makes detailed output to stdout
 * `--silent, -s` - makes no output to stdout. Overrides `-v` argument

Typical blockfinder output:
~~~~
[NCS = NC2]
[ELB samples = 2 patterns = 2]
NC
CN
[ELB samples = 3 patterns = 3]
NNC
NCN
CNN
[ELB samples = 4 patterns = 4]
NNNC
NNCN
NCNN
CNNN
...
~~~~
---
### elbstats.py

**elbstats.py** takes ELB-file as input and calculates how many blocks
of each size are contained in the given ELB-file.

Typical **elbstats.py** output:
~~~~
Blocks used:
 2 X  2 - 1 block(s)
 3 X  3 - 1 block(s)
 4 X  4 - 2 block(s)
 5 X  5 - 4 block(s)
 5 X  6 - 1 block(s)
 6 X  7 - 10 block(s)
 6 X  8 - 2 block(s)
 6 X  9 - 1 block(s)
 7 X 10 - 8 block(s)
 7 X 11 - 1 block(s)
-----------
 2 X . - 1 block(s)
 3 X . - 1 block(s)
 4 X . - 2 block(s)
 5 X . - 5 block(s)
 6 X . - 13 block(s)
 7 X . - 9 block(s)
 TOTAL - 31 blocks
~~~~
----
### elbclean.py

**elbclean.py** takes ELB-file as input and cleans redundant blocks from
the file. Writes separate "clean" ELB-file.

---
### UUCSL.py

**UUCSL.py** is the main script in the package. It calculates the
resulting price-optimized labeling schemes using direct products of
ELBs (concatenated patterns of two or more ELBs in all possible
combinations) and simplex algorithm.

It has 3 positional command-line arguments:
* `ncs` - NCS
* `elb_file` - ELB-file
* `price_file` - file with prices (prices of commercially available in
2016 amino acids with various labeling types are collected in `csv`
prices-files in `prices` folder of the package)

Optional command line arguments:
* `--aminoacids, -a` - list of amino acids (one-letter code) in the
 final scheme. By default the script uses the list of all 20 biogenic
 amino acids.
* `--notaminoacids, -n` - list of amino acids that will not be
 concidered in scheme calculation. This command will override
 intersections with `-a` argument.
* `--jobname, -j` - this value will be used as a prefix to all output files
* `--verbose, -v` - makes detailed output to stdout
* `--silent, -s` - makes no output to stdout. Overrides `-v` argument

Typical output of **UUCSL.py**

~~~~
Best scheme
[NCS = NCD4]
[solution]
Res,S1,S2,S3,S4,S5,S6
Ala, N, D, D, N, N, C
Arg, D, X, D, N, D, D
Asn, D, X, D, D, N, D
Asp, N, D, D, N, D, D
Cys, X, N, N, N, N, C
Gln, D, X, D, N, N, C
Glu, N, D, D, D, N, D
Gly, N, D, D, N, C, N
His, X, N, N, N, D, D
Ile, D, X, D, N, C, N
Leu, D, X, D, C, N, N
Lys, D, D, X, N, D, D
Met, X, N, N, N, C, N
Phe, D, D, X, N, N, C
Pro, X, X, X, F, X, F
Ser, D, D, X, N, C, N
Thr, D, D, X, D, N, D
Trp, X, N, N, C, N, N
Tyr, D, D, X, C, N, N
Val, N, D, D, C, N, N
[price]
13530.304579999996

Blocks used for this scheme:
[ELB samples = 3 patterns = 4]
XNN
NDD
DXD
DDX
[ELB samples = 3 patterns = 5]
NNC
NCN
NDD
CNN
DND
~~~~
---
### CombLabel_check.py

**CombLabel_check.py** - calculates labeling dictionary and assignment
statistics given NCS, labeling scheme and protein sequence.

Positional arguments:
* `ncs` - NCS (NMR coding system)
* `sequence` - protein sequence (one-letter code, `*.fasta` files)
* `solution` - solution file, e.g. output file of `UUCSL.py`

Optional argument:
* `--prices, -p` - prices-file, like one that can be used by `UUCSL.py`

Example of labeling dictionary (`*_dictionary.txt`):
~~~~
...
0112110: P60 - Q61
0121021: L8 - W9
0121120: L120 - Q121
0121120: L150 - Q151
0121120: L47 - Q48
0121210: F122 - Q123
0122011: K87 - W88
0211102: I134 - M135
0211102: I29 - M30
0211210: V129 - Q130
0212201: T74 - M75
0221021: A71 - W72
...
~~~~

Assignment statistics (`*_solution.txt`) contains a variety of
statistical information about protein sequence and given solution.

---

## Example of usage

to be added


## Future improvements:

* make `blockfinder.py` and `UUCSL.py` calculations parallel
* code refactoring
* precompilation of core components to improve speed

## Citation

When using this package for scientific use, please cite the following
paper:

M.Dubinnyi, M.Myshkin, Z.Shenkarev: "UUCSL: Universal unambiguous
combinatorial selective labeling schemes" JBNMR, 2019