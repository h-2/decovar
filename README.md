# deCoVar – the 𝐝𝐞𝐂𝐎DE 𝐯𝐚𝐫iant tools

A application with multiple sub-programs that perform useful transformations, filters, etc. on VCF/BCF files.
Developed at deCODE genetics to fill-in gaps of bcftools.

## Sub-programs

* `allele_shrink`: reduce the size impact of multiallelic records by removing rare alleles and/or replacing the `PL`
field with

## Notable differences to BCFtools

* No 2GB limit for records when reading/writing `.vcf` and `.vcf.gz` (we can do nothing about `.bcf` unfortunately).

## Instructions

### Getting it

We recommed using a package from the [releases page](./releases).

<details>
<summary>Cloning and building from source</summary>

<p>

*Please note that GCC>=10.3 is required; LLVM and MSVC are not supported.*

Clone the latest main-branch:

```
cd ~/devel/                                                         # or some other directory
git clone --recurse-submodules https://github.com/h-2/decovar.git
```

Setup build folder:

```
mkdir -p ~/devel/decovar-build/release                              # or some other directory
cd ~/devel/decovar-build/release
cmake -DCMAKE_BUILD_TYPE=Release ../../decovar
```

Build:

```
make -j 4
```

Run:

```
./bin/decovar --help
```
</p>
</details>

### Using

Look at the help-page for details:

```
path/to/decovar --help
```

## Disclaimer

* This is an early preview and everything is still subject to change.
