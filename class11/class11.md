# Class 11: AlphaFold
Melissa (PID: A19149673)

- [Custom Analysis of Resulting
  Models](#custom-analysis-of-resulting-models)
- [Predicted alignment error for
  domains](#predicted-alignment-error-for-domains)
- [Residue conservation from alignment
  file](#residue-conservation-from-alignment-file)

## Custom Analysis of Resulting Models

``` r
library(bio3d)
pdb <- read.pdb("hivp_dimer_0d832/hivp_dimer_0d832_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb")
pdb
```


     Call:  read.pdb(file = "hivp_dimer_0d832/hivp_dimer_0d832_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb")

       Total Models#: 1
         Total Atoms#: 1440,  XYZs#: 4320  Chains#: 1  (values: A)

         Protein Atoms#: 1440  (residues/Calpha atoms#: 176)
         Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)

         Non-protein/nucleic Atoms#: 0  (residues: 0)
         Non-protein/nucleic resid values: [ none ]

       Protein sequence:
          MRDEIATTVFFVTRLVKKHNKLSKQIEDFAEKLMTILFETYRSHWHSDCPSKGQAFRCIR
          INNNQNKDPIERACAESKVDFFHLGLPKELTIWDPFEVCCRYGEKSHPFTVASFKGRWEE
          WELYQQISYAVSRASSDDSSGTSCNEESCGKEPHVIPKVSNPKSIYQVENLKQPFQ

    + attr: atom, xyz, calpha, call

Make a vector of input PDB file names that we can read into R

``` r
pdb_files <- list.files("hivp_dimer_0d832/", 
           pattern =".pdb",
           full.names = TRUE)
basename(pdb_files)
```

    [1] "hivp_dimer_0d832_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
    [2] "hivp_dimer_0d832_unrelaxed_rank_002_alphafold2_ptm_model_2_seed_000.pdb"
    [3] "hivp_dimer_0d832_unrelaxed_rank_003_alphafold2_ptm_model_3_seed_000.pdb"
    [4] "hivp_dimer_0d832_unrelaxed_rank_004_alphafold2_ptm_model_4_seed_000.pdb"
    [5] "hivp_dimer_0d832_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000.pdb"

``` r
# Read all data from Models 
#  and superpose/fit coords
pdbs <- pdbaln(pdb_files, fit=TRUE, exefile="msa")
```

    Reading PDB files:
    hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
    hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_002_alphafold2_ptm_model_2_seed_000.pdb
    hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_003_alphafold2_ptm_model_3_seed_000.pdb
    hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_004_alphafold2_ptm_model_4_seed_000.pdb
    hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000.pdb
    .....

    Extracting sequences

    pdb/seq: 1   name: hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb 
    pdb/seq: 2   name: hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_002_alphafold2_ptm_model_2_seed_000.pdb 
    pdb/seq: 3   name: hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_003_alphafold2_ptm_model_3_seed_000.pdb 
    pdb/seq: 4   name: hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_004_alphafold2_ptm_model_4_seed_000.pdb 
    pdb/seq: 5   name: hivp_dimer_0d832//hivp_dimer_0d832_unrelaxed_rank_005_alphafold2_ptm_model_5_seed_000.pdb 

``` r
pdbs
```

                                   1        .         .         .         .         50 
    [Truncated_Name:1]hivp_dimer   MRDEIATTVFFVTRLVKKHNKLSKQIEDFAEKLMTILFETYRSHWHSDCP
    [Truncated_Name:2]hivp_dimer   MRDEIATTVFFVTRLVKKHNKLSKQIEDFAEKLMTILFETYRSHWHSDCP
    [Truncated_Name:3]hivp_dimer   MRDEIATTVFFVTRLVKKHNKLSKQIEDFAEKLMTILFETYRSHWHSDCP
    [Truncated_Name:4]hivp_dimer   MRDEIATTVFFVTRLVKKHNKLSKQIEDFAEKLMTILFETYRSHWHSDCP
    [Truncated_Name:5]hivp_dimer   MRDEIATTVFFVTRLVKKHNKLSKQIEDFAEKLMTILFETYRSHWHSDCP
                                   ************************************************** 
                                   1        .         .         .         .         50 

                                  51        .         .         .         .         100 
    [Truncated_Name:1]hivp_dimer   SKGQAFRCIRINNNQNKDPIERACAESKVDFFHLGLPKELTIWDPFEVCC
    [Truncated_Name:2]hivp_dimer   SKGQAFRCIRINNNQNKDPIERACAESKVDFFHLGLPKELTIWDPFEVCC
    [Truncated_Name:3]hivp_dimer   SKGQAFRCIRINNNQNKDPIERACAESKVDFFHLGLPKELTIWDPFEVCC
    [Truncated_Name:4]hivp_dimer   SKGQAFRCIRINNNQNKDPIERACAESKVDFFHLGLPKELTIWDPFEVCC
    [Truncated_Name:5]hivp_dimer   SKGQAFRCIRINNNQNKDPIERACAESKVDFFHLGLPKELTIWDPFEVCC
                                   ************************************************** 
                                  51        .         .         .         .         100 

                                 101        .         .         .         .         150 
    [Truncated_Name:1]hivp_dimer   RYGEKSHPFTVASFKGRWEEWELYQQISYAVSRASSDDSSGTSCNEESCG
    [Truncated_Name:2]hivp_dimer   RYGEKSHPFTVASFKGRWEEWELYQQISYAVSRASSDDSSGTSCNEESCG
    [Truncated_Name:3]hivp_dimer   RYGEKSHPFTVASFKGRWEEWELYQQISYAVSRASSDDSSGTSCNEESCG
    [Truncated_Name:4]hivp_dimer   RYGEKSHPFTVASFKGRWEEWELYQQISYAVSRASSDDSSGTSCNEESCG
    [Truncated_Name:5]hivp_dimer   RYGEKSHPFTVASFKGRWEEWELYQQISYAVSRASSDDSSGTSCNEESCG
                                   ************************************************** 
                                 101        .         .         .         .         150 

                                 151        .         .     176 
    [Truncated_Name:1]hivp_dimer   KEPHVIPKVSNPKSIYQVENLKQPFQ
    [Truncated_Name:2]hivp_dimer   KEPHVIPKVSNPKSIYQVENLKQPFQ
    [Truncated_Name:3]hivp_dimer   KEPHVIPKVSNPKSIYQVENLKQPFQ
    [Truncated_Name:4]hivp_dimer   KEPHVIPKVSNPKSIYQVENLKQPFQ
    [Truncated_Name:5]hivp_dimer   KEPHVIPKVSNPKSIYQVENLKQPFQ
                                   ************************** 
                                 151        .         .     176 

    Call:
      pdbaln(files = pdb_files, fit = TRUE, exefile = "msa")

    Class:
      pdbs, fasta

    Alignment dimensions:
      5 sequence rows; 176 position columns (176 non-gap, 0 gap) 

    + attr: xyz, resno, b, chain, id, ali, resid, sse, call

``` r
rd <- rmsd(pdbs, fit=T)
```

    Warning in rmsd(pdbs, fit = T): No indices provided, using the 176 non NA positions

``` r
range(rd)
```

    [1]  0.000 21.947

``` r
library(pheatmap)

colnames(rd) <- paste0("m",1:5)
rownames(rd) <- paste0("m",1:5)
pheatmap(rd)
```

![](Class11_files/figure-commonmark/unnamed-chunk-5-1.png)

The RMSD heatmap compares structural similarity between the five
AlphaFold models. Additionally, the lower RMSD values indicate more
similar structures while models with higher RMSD are more structurally
different from the others.

``` r
pdb <- read.pdb("1hsg")
```

      Note: Accessing on-line PDB file

``` r
plotb3(pdbs$b[1,], typ="l", lwd=2, sse=pdb)
```

    Warning in plotb3(pdbs$b[1, ], typ = "l", lwd = 2, sse = pdb): Length of input
    'sse' does not equal the length of input 'x'; Ignoring 'sse'

``` r
points(pdbs$b[2,], typ="l", col="red")
points(pdbs$b[3,], typ="l", col="blue")
points(pdbs$b[4,], typ="l", col="darkgreen")
points(pdbs$b[5,], typ="l", col="orange")
abline(v=100, col="gray")
```

![](Class11_files/figure-commonmark/unnamed-chunk-7-1.png)

The pLDDT plot shows AlphaFold confidence across residue positions.
Moreover, the higher pLDDT values indicate higher confidence, while
lower values suggest less reliable regions of the model.

``` r
core <- core.find(pdbs)
```

     core size 175 of 176  vol = 59740.11 
     core size 174 of 176  vol = 49331.57 
     core size 173 of 176  vol = 46442.46 
     core size 172 of 176  vol = 44034.47 
     core size 171 of 176  vol = 41663.07 
     core size 170 of 176  vol = 39972.2 
     core size 169 of 176  vol = 37971.21 
     core size 168 of 176  vol = 36135.77 
     core size 167 of 176  vol = 34793.5 
     core size 166 of 176  vol = 33767.56 
     core size 165 of 176  vol = 32341.92 
     core size 164 of 176  vol = 31045.85 
     core size 163 of 176  vol = 29234.63 
     core size 162 of 176  vol = 27439.82 
     core size 161 of 176  vol = 25645.37 
     core size 160 of 176  vol = 23902.96 
     core size 159 of 176  vol = 21973.68 
     core size 158 of 176  vol = 20161.16 
     core size 157 of 176  vol = 18486.13 
     core size 156 of 176  vol = 17019.19 
     core size 155 of 176  vol = 15769.75 
     core size 154 of 176  vol = 14503.63 
     core size 153 of 176  vol = 13657.27 
     core size 152 of 176  vol = 13029.01 
     core size 151 of 176  vol = 12707.01 
     core size 150 of 176  vol = 12781.06 
     core size 149 of 176  vol = 12943.19 
     core size 148 of 176  vol = 12555.66 
     core size 147 of 176  vol = 12205.02 
     core size 146 of 176  vol = 12247.06 
     core size 145 of 176  vol = 11547.95 
     core size 144 of 176  vol = 10764.27 
     core size 143 of 176  vol = 9999.411 
     core size 142 of 176  vol = 9230.039 
     core size 141 of 176  vol = 8041.899 
     core size 140 of 176  vol = 6786.012 
     core size 139 of 176  vol = 5173.183 
     core size 138 of 176  vol = 3247.493 
     core size 137 of 176  vol = 911.866 
     core size 136 of 176  vol = 210.723 
     core size 135 of 176  vol = 162.178 
     core size 134 of 176  vol = 138.219 
     core size 133 of 176  vol = 118.938 
     core size 132 of 176  vol = 106.287 
     core size 131 of 176  vol = 95.029 
     core size 130 of 176  vol = 81.78 
     core size 129 of 176  vol = 73.211 
     core size 128 of 176  vol = 64.041 
     core size 127 of 176  vol = 56.801 
     core size 126 of 176  vol = 47.741 
     core size 125 of 176  vol = 39.084 
     core size 124 of 176  vol = 33.14 
     core size 123 of 176  vol = 28.014 
     core size 122 of 176  vol = 22.268 
     core size 121 of 176  vol = 17.934 
     core size 120 of 176  vol = 13.563 
     core size 119 of 176  vol = 9.068 
     core size 118 of 176  vol = 4.391 
     core size 117 of 176  vol = 1.548 
     core size 116 of 176  vol = 0.652 
     core size 115 of 176  vol = 0.523 
     core size 114 of 176  vol = 0.355 
     FINISHED: Min vol ( 0.5 ) reached

``` r
core.inds <- print(core, vol=0.5)
```

    # 115 positions (cumulative volume <= 0.5 Angstrom^3) 
      start end length
    1     1  65     65
    2    67 116     50

``` r
xyz <- pdbfit(pdbs, core.inds, outpath="corefit_structures")
```

``` r
rf <- rmsf(xyz)

plotb3(rf, sse=pdb)
```

    Warning in plotb3(rf, sse = pdb): Length of input 'sse' does not equal the
    length of input 'x'; Ignoring 'sse'

``` r
abline(v=100, col="gray", ylab="RMSF")
```

![](Class11_files/figure-commonmark/unnamed-chunk-11-1.png)

The RMSF plot shows how much each residue varies across the fitted
models as the higher RMSF values indicate more flexible or variable
regions.

## Predicted alignment error for domains

``` r
library(jsonlite)
# Listing of all PAE JSON files
results_dir <- "hivp_dimer_0d832/"
pae_files <- list.files(path=results_dir,
                        pattern=".*model.*\\.json",
                        full.names = TRUE)
```

``` r
pae1 <- read_json(pae_files[1],simplifyVector = TRUE)
pae5 <- read_json(pae_files[5],simplifyVector = TRUE)

attributes(pae1)
```

    $names
    [1] "plddt"   "max_pae" "pae"     "ptm"    

``` r
# Per-residue pLDDT scores 
#  same as B-factor of PDB..
head(pae1$plddt)
```

    [1] 87.50 90.62 91.25 92.25 94.38 95.44

``` r
pae1$max_pae
```

    [1] 31.29688

``` r
pae5$max_pae
```

    [1] 31.14062

The maximum PAE values compare prediction uncertainty between models. A
lower max PAE indicates a more reliable model. Hence, the model 1 is
more reliable than the model 5.

``` r
plot.dmat(pae1$pae, 
          xlab="Residue Position (i)",
          ylab="Residue Position (j)")
```

![](Class11_files/figure-commonmark/unnamed-chunk-16-1.png)

``` r
plot.dmat(pae5$pae, 
          xlab="Residue Position (i)",
          ylab="Residue Position (j)",
          grid.col = "black",
          zlim=c(0,30))
```

![](Class11_files/figure-commonmark/unnamed-chunk-17-1.png)

``` r
plot.dmat(pae1$pae, 
          xlab="Residue Position (i)",
          ylab="Residue Position (j)",
          grid.col = "black",
          zlim=c(0,30))
```

![](Class11_files/figure-commonmark/unnamed-chunk-18-1.png)

The PAE plots show predicted alignment error between residue pairs. The
first Model 1 shows lower overall error than model 5. This suggests a
more confident domain arrangement.

## Residue conservation from alignment file

``` r
aln_file <- list.files(path=results_dir,
                       pattern=".a3m$",
                        full.names = TRUE)
aln_file
```

    [1] "hivp_dimer_0d832//hivp_dimer_0d832.a3m"

``` r
aln <- read.fasta(aln_file[1], to.upper = TRUE)
```

    [1] " ** Duplicated sequence id's: 101 **"

``` r
dim(aln$ali)
```

    [1] 2704  232

``` r
sim <- conserv(aln)
plotb3(sim[1:99], sse=trim.pdb(pdb, chain="A"),
       ylab="Conservation Score")
```

![](Class11_files/figure-commonmark/unnamed-chunk-22-1.png)

The conservation plot shows which residues are most conserved across the
sequence alignment.

``` r
con <- consensus(aln, cutoff = 0.9)
con$seq
```

      [1] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [19] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [37] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [55] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [73] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [91] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [109] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [127] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [145] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [163] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [181] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [199] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [217] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"

``` r
m1.pdb <- read.pdb(pdb_files[1])

# uses conservation scores for all residues in the model
sim2 <- rep(sim[1:99], 2)

# map residue conservation scores onto every atom
occ <- sim2[match(m1.pdb$atom$resno, sort(unique(m1.pdb$atom$resno)))]

# write new PDB file for Mol Star
write.pdb(m1.pdb, o = occ, file = "m1_conserv.pdb")
```

The original code produced a length mismatch error, so it was adjusted
to map conservation scores to residue numbers.
