Archaic introgression lecture exercises (February 2025)
================

## Exercise 1 – testing introgression using the $f_4$ statistic

If you want to see this document with all the solutions included, click
[here](README_solutions.md) (but please try to solve the exercises on
your own at first!).

You can find the slides from the lecture
[here](https://github.com/bodkan/ku-introgression2025/blob/main/lecture.pdf).

### Introduction

*You sequenced the genomes of four Africans and four Eurasians and got
genotypes from a single chromosome from each of them (i.e., you have the
genotypes of four African and four Eurasian chromosomes in total).
Unfortunately, there’s been a mix up in the lab and you don’t know which
one of them is African and which is Eurasian! You only know that they
are labeled A, B, C, …, H. What a disaster!*

*Fortunately, you also have genotypes from three other individuals whose
identity you know for certain: an African, a Neanderthal, and a
chimpanzee. This means that you are able to compute an* $f_4$ statistic
which will test for evidence of Neanderthal introgression in a given
some sample $X$.

*Can you save the day and determine which of the A, B, C, …, H samples
are African and which are Eurasian based on the following* $f_4$
statistic test?

$$
f_4(\textrm{African}, X; \textrm{Neanderthal}, \textrm{Chimp}).
$$

*Recall that only Eurasians are expected to have appreciable amounts of
Neanderthal ancestry but Africans don’t. So, hopefully, by computing
this statistic by swapping for X all of the mixed up samples A, B, C…,
you should be able to determine which are African (these won’t show
evidence of Neanderthal introgression, giving* $f_4$ *values close to
zero)* *and which are Eurasian (these will give* $f_4$ *values
significantly more negative).*

### Moving over to R

Type `R` in your terminal or (better) just use RStudio R console if you
have it.

#### Task: Read and inspect the genotypes of all the sequenced samples

First **read the genotype table into R**:

``` r
gt <- read.table(url("https://github.com/bodkan/ku-introgression2025/raw/main/genotypes_ex1.tsv"), sep = "\t", header = TRUE)
```

**Familiarize yourself with the data** by running this R command which
shows genotype information from only the first few sites in the genome:

``` r
head(gt)
```

The `gt` data set is a plain R data frame where each column contains the
genotype of that individual (`0` - ancestral allele, `1` - derived
allele).

#### Task: Count SNPs

**For how many loci in the genome do we have genotype data available?**

``` r
nrow(gt)
```

#### Task: Count AFR-Chimp, NEA-Chimp, AFR-NEA shared alleles

You can extract all genotypes of a given individual as a single vector
by using the `$` or `[[` subsetting operators of R data frames like
this:

``` r
gt$African

gt[["African"]]
```

A useful trick for comparing two chromosomes in their entirety is to
rely on the fact that R can perform *vectorized operations* (operations
performed on multiple elements of a vector at once). For instance, if
this gives us the genotypes of an African and Neanderthal chromosome:

``` r
gt[["African"]]

gt[["Neanderthal"]]
```

then we can find positions at which those two samples carry the same
allele like this:

``` r
gt[["African"]] == gt[["Neanderthal"]] # this gives us TRUE/FALSE values 
```

Counting those positions can be done using the `sum()` function like
this:

``` r
# sum() treats TRUE as 1 and FALSE as 0, so we can sum everything up!
# -- this gives us the number of positions at which an African carries the same
#    allele as the Neanderthal
sum(gt[["African"]] == gt[["Neanderthal"]])
```

So the answer to this task’s question can be computed as:

``` r
sum(gt[["African"]] == gt[["Chimp"]])
sum(gt[["Neanderthal"]] == gt[["Chimp"]])
sum(gt[["African"]] == gt[["Neanderthal"]])
```

**Do the counts of allele sharing that you got make sense from a
phylogenetic point of view?**

#### Task: Compute $f_4(\textrm{AFR, X; NEA, Chimp})$ for one of the unknown samples A-H

Above we computed alleles which *agree* between two samples.

On the other hand, this would count how many alleles are *different*
between a African and chimpanzee chromosome:

``` r
sum(gt[["African"]] != gt[["Chimp"]]) # note the != instead of ==
```

Inside the `sum()` function we can compose multiple logical conditions
to create more complex comparison operations using the `&` operator (AND
operation in mathematical logic).

Armed with this knowledge, we can compute the BABA and ABBA counts using
this bit of R code:

``` r
X <- "A"

abba <- sum(
  (gt[["African"]] == gt[["Chimp"]]) &         # filters for A**A sites
  (gt[["African"]] != gt[["Neanderthal"]]) &   # filters for A*B* sites
  (gt[[X]]         == gt[["Neanderthal"]])     # filters for *BB* sites
)                                              # together then ABBA

baba <- sum(
  (gt[["African"]] != gt[["Chimp"]]) &         # filters for B**A sites
  (gt[["African"]] == gt[["Neanderthal"]]) &   # filters for B*B* sites
  (gt[[X]]         == gt[["Chimp"]])           # filters for *A*A sites
)                                              # together then BABA
```

From these counts we can get an idea about whether one or the other are
more frequently appearing in the data:

``` r
baba - abba
```

Finally, we can compute an $f_4$ statistic like this, which simply
normalizes the raw difference between BABA and ABBA counts by how many
SNPs we have in our data set:

``` r
f4_value <- (baba - abba) / nrow(gt)

f4_value
```

#### Task: Are ABBA or BABA sites the only ones in our data? For instance, can you find if there are any AAAB sites for the quartet $f_4(\textrm{AFR, X; NEA, Chimp})$? Would those be useful for studying introgression?

``` r
X <- "A"

aaab <- sum(
  (gt[["African"]]     == gt[[X]]) &               # filters for AA** sites
  (gt[[X]]             == gt[["Neanderthal"]]) &   # filters for *AA* sites
  (gt[["Neanderthal"]] != gt[["Chimp"]])           # filters for **AB sites
)

aaab
```

#### Task (full solution under the line below):

**You know that if `X` is a African, you expect to see roughly the same
count of `BABA` and `ABBA` site patterns from the**
$f_4(\textrm{AFR, X; NEA, Chimp})$ **quartet, so the difference should
“be about zero”. Use the code above to compute `baba`, `abba`, and
`f4_value` to for all of your mixed up samples A, B, C, …, H and note
down the values you got for each – which samples are most likely African
and which ones are Eurasian?**

*\[If you are more familiar with R, compute the counts automatically in
a loop of some kind and make a figure.\]*

#### Task: **What does it mean for this test statistic to “be about zero”? What are we missing to truly use this as a statistical significance test?**

------------------------------------------------------------------------

If you’re not comfortable with R, feel free to run this in full and
answer the questions from tasks above based on the results you get:

``` r
X <- c("A", "B", "C", "D", "E", "F", "G", "H")

f4_values <- sapply(X, function(x) {
  abba <- sum(
    (gt[["African"]] == gt[["Chimp"]]) &         # filters for A**A sites
    (gt[["African"]] != gt[["Neanderthal"]]) &   # filters for A*B* sites
    (gt[[x]]         == gt[["Neanderthal"]])     # filters for *BB* sites
  )                                              # together then ABBA
  
  baba <- sum(
    (gt[["African"]] != gt[["Chimp"]]) &         # filters for B**A sites
    (gt[["African"]] == gt[["Neanderthal"]]) &   # filters for B*B* sites
    (gt[[x]]         == gt[["Chimp"]])           # filters for *A*A sites
  )                                              # together then BABA
  
  (baba - abba) / nrow(gt)
})

data.frame(f4_values) # to print values in a neater table
```

``` r
plot(f4_values, xaxt = "n", xlab = "test sample", ylab = "f4(African, X; Neanderthal, Chimp)")
abline(h = 0, lty = 2, col = "red")
axis(side = 1, at = seq_along(X), labels = X)
```

We can see that the samples A-D are consistent with an $f_4$ statistic
“value of about 0”“, meaning that the BABA and ABBA counts were”about
the same”. This is what we would expect for African samples who are not
expected to be closer to a Neanderthal genome than another African.

On the other hand, samples E-H show a much “more negative value of the
$f_4$ statistic”“, which is consistent with an excess of ABBA sites
compared to BABA sites – which arise with an increased sharing of
derived alleles between the sample X and a Neanderthal genome, just as
we would expect when X is of Eurasian ancestry.

**Important:** In this simple example we’re missing confidence intervals
– those would allow us to do a proper statistical test to determine for
which samples we really cannot reject a null hypothesis of no gene flow
from Neanderthals. Having the information about confidence intervals
would, we could avoid the vague and statistically unsatisfying talk
about some value being “almost zero”, and some other value being “much
more negative” than that (in this exercise we did this for simplicity).
The confidence interval for a given $f_4$ statistic would either
intersect the 0 null hypothesis or not. For an example, see Figure 3 in
[this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6485383/#F3).

Real-world software such as
[ADMIXTOOLS](https://github.com/DReichLab/AdmixTools) computes
confidence intervals using a so-called
[bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
procedure across windows along a genome.

------------------------------------------------------------------------

If you want to take a closer look at how the genotype data was prepared
(it was simulated!), you can see the complete code
[here](generate_genotypes.R).

## Exercise 2 – estimating the proportion of Neanderthal ancestry

### Introduction

*Having saved the day by identifying which of the A-H samples are of
African or Eurasian origin by testing which of them appear to carry
evidence of Neanderthal introgression, you now want to estimate how much
of their genome derives from the Neanderthals. In order to do this, you
need to compute a ratio of* $f_4$ *values as described in the lecture.*

*Of course, in order to compute this* $f_4$*-ratio estimate, you will
need “another Neanderthal” genome! Luckily, we now have genomes of
several Neanderthals so this is not an issue and a local friendly
bioinformatician has already presciently merged your `gt` genotype table
from the first exercise with the genotypes of “another Neanderthal”.*

*Estimate the proportion of Neanderthal ancestry in each of your A-H
samples!*

### Moving over to R

Type `R` in your terminal or (better) just use RStudio R console if you
have it.

#### Task: Read and inspect the genotypes of all the sequenced samples

**You will be using the same genotype table as in the previous exercise,
with one additional column called `another_Neanderthal`. You can read it
again like this:**

``` r
gt <- read.table(url("https://github.com/bodkan/ku-introgression2025/raw/main/genotypes_ex2.tsv"), sep = "\t", header = TRUE)
```

As always, **verify that the format of the data set you have matches
what you expect:**

``` r
head(gt)

nrow(gt)
```

#### Task: Estimate Neanderthal ancestry proportion in samples A-H

From the lecture you know that we can get an estimate for the proportion
of Neanderthal ancestry in a sample $X$ by dividing the rate of allele
sharing between $X$ and a Neanderthal genome (one $f_4$ statistic) by
the rate of allele sharing expected between two Neanderthals (another
$f_4$ statistic).

To do this, we can take the $f_4$ values you computed in Exercise 1 for
all samples A-H, and divide those values by
$f_4(\textrm{African, another Neanderthal; Neanderthal, Chimp})$ (which
we can do with the new set of genotypes `gt`):

``` r
# we can compute the f4 values for everyone (A-H samples as well as
# "another_neanderthal") using the same code as above
X <- c("A", "B", "C", "D", "E", "F", "G", "H",
       "another_Neanderthal") # <---- we added this to our loop

f4_values <- sapply(X, function(x) {
  abba <- sum(
    (gt[["African"]] == gt[["Chimp"]]) &         # filters for A**A sites
    (gt[["African"]] != gt[["Neanderthal"]]) &   # filters for A*B* sites
    (gt[[x]]         == gt[["Neanderthal"]])     # filters for *BB* sites
  )                                              # together then ABBA
  
  baba <- sum(
    (gt[["African"]] != gt[["Chimp"]]) &         # filters for B**A sites
    (gt[["African"]] == gt[["Neanderthal"]]) &   # filters for B*B* sites
    (gt[[x]]         == gt[["Chimp"]])           # filters for *A*A sites
  )                                              # together then BABA
  
  (baba - abba) / nrow(gt)
})


# to arrive at the estimate of Neanderthal ancestry, we divide f4 values for
# samples A-H by f4 value comparing the two Neanderthals
proportions <- f4_values / f4_values["another_Neanderthal"]

data.frame(proportions) # to print values in a neater table
```

#### Task: Plot the estimated proportions of Neanderthal ancestry

**How much Neanderthal ancestry did you estimate in Africans vs
Eurasians? Do those numbers fit what you’ve learned from the lecture?**

To make the results clearer to see, let’s visualize them:

``` r
plot(proportions[-length(proportions)] * 100,
     xlab = "test individual", ylab = "proportion of Neanderthal ancestry [%]",
     ylim = c(0, 10))

abline(h = 3, lty = 2, col = "red")
```

#### Task: **Why did we not plot the proportion of Neanderthal ancestry in the very last item of the `proportions` variable? What does that last element of the vector `proportions` contain and why?**

------------------------------------------------------------------------

If you want to take a closer look at how the genotype data was prepared
(it was simulated!), you can see the complete code
[here](generate_genotypes.R).

<!--
## Exercise 3 -- dating Neanderthal admixture
&#10;(This exercise is a bonus for those of you who are already experts in R. If you're not very comfortable with R or population genetics, don't feel like you have to rush through the exercises to finish everything. You can take a look at this exercise at home, if you want.)
&#10;### Introduction
&#10;*You sequenced 100 diploid genomes from a Eurasian population and are interested in estimating the time of Neanderthal introgression into the ancestors of this population. The literature suggests that the introgression happened around 55 thousand years ago -- does this also apply to the population that you sequenced?*
&#10;*To be able to do this, you ran an inference software which gives you the exact coordinates of Neanderthal DNA tracts present in every Eurasian genome that you sequenced. This of course means that you also know the lengths of each of those tracts.*
&#10;*Use the distribution of the Neanderthal tract lengths in your population to estimate the time of Neanderthal introgression!*
&#10;### Moving over to R
&#10;Type `R` in your terminal or (better) just use RStudio R console if you have it.
&#10;#### Task: Load and inspect the tracts data
&#10;**First load the table with coordinates of all Neanderthal tracts** into R:
&#10;
``` r
tracts <- read.table(url("https://github.com/bodkan/ku-introgression2025/raw/main/tracts.tsv"), sep = "\t", header = TRUE)
```
&#10;**Familiarize yourself with the data** by running this R command which shows information from only the first few genotypes:
&#10;
``` r
head(tracts)
```
&#10;**For how many individuals do we have information about Neanderthal DNA tracts that they carry?**
&#10;
``` r
length(unique(tracts$individual))
```
&#10;#### Task: plot the distribution of tract lengths across bins
&#10;It looks like the inference software (or a helpful bioinformatician) binned each tract according to its length (see the column `bin`). **What does the distribution of Neanderthal tract lengths looks like in your data?** Knowing that recombination has acted on the introgressed Neanderthal DNA over time, each generation, suggests that the distribution should look exponential -- do you see this in the data? To answer this, plot the proportion of tracts in each bin.
&#10;
``` r
# get the bin numbers
bins <- sort(unique(tracts$bin))
&#10;# count the tracts in each bin and compute the proportion of tracts in each bin
counts <- as.integer(table(tracts$bin))
props <- counts / sum(counts)
&#10;plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion of tracts")
```
&#10;The distribution does, indeed, look quite exponential. **Next we will try to use this to date Neanderthal introgression using the information encoded in the distribution of tract lengths!**
&#10;#### Dating the introgression event -- a bit of theory first
&#10;As we know, over time since admixture, recombination breaks up longer haplotypes into shorter ones, regularly almost like a clock. And it turns out that the distribution of tract lengths after time $t$ follows exponential decay, leading to the distribution of tract lengths $x$ to have the following form:
&#10;$$
\textrm{tract length of }~x \sim \exp(-\lambda x) = \exp(-r t x),
$$
&#10;where the $\lambda$ parameter determines the *rate* of exponential decay which can be computed, under some simplifying assumptions, as the product of the recombination rate (traditionally in humans with value of about $1 \times 10^{-8}$) and $t$ which is the time since admixture:
&#10;$$
\lambda = rt
$$
&#10;**The** $t$ in this equation latter is our unknown we're trying to compute in this exercise! So we know which equation we can use to extimate the admixture time.
&#10;It also turns out that the expected value of this exponential distribution plotted above (which [can be computed](https://en.wikipedia.org/wiki/Exponential_distribution#Mean,_variance,_moments,_and_median) simply as $1 / \lambda = 1 / rt$) gives us the theoretical expression for the expected tract length after time $t$.
&#10;#### Task: Compute the average length of an introgressed tract
&#10;Of course, you can also compute this expectation from the data, **so do this now: get an estimate of the expected length of an introgressed fragment after (unknown) time** $t$ by computing the average introgressed tract length observed in data. What is the average length of a Neanderthal DNA segment in the data? (Remember that at the moment of introgression, *entire Neanderthal chromosomes* were segregating in modern humans!)
&#10;
``` r
L <- mean(tracts$length)
L # length in units of base pairs
```
&#10;#### Task: date the introgression event -- in practice
&#10;Taking the simple math above together and doing a little algebra, we can express the average expected length of an introgressed fragment after time $t$ using this formula:
&#10;$$
\textrm{average tract length}~L = \frac{1}{\lambda} = \frac{1}{rt}
$$
&#10;**Because we know** $L$ as computed just above (average tract length as you just computed) and $r$ (recombination rate of $1 \times 10^{-8}$), this means we can estimate the time since the admixture by a simple rearrangement of the above equation to separate the time of introgression $t$:
&#10;$$
t = \frac{1}{rL},
$$
&#10;where `r` is the recombination rate and `L` is the average length of an introgressed Neanderthal tract (as you just computed it).
&#10;Note that this time estimate will be in units of generations, so we'll have to multiply this quantity by the generation time (roughly 30 years for humans) to get time in years before the present.
&#10;**Use this simple equation** $\frac{1}{rL}$ to compute the estimate of the admixture time:
&#10;
``` r
r <- 1e-8 # crossovers per bp per generation
L <- mean(tracts$length) # average tract length after time t
&#10;t <- 1 / (r * L)
t
&#10;# convert the time of introgression to years before present assuming
# generation time of 30 years
t * 30
```
&#10;You should get a value will be quite close to \~55 thousand years ago, an estimate which is often found in the literature as the time when Neanderthals and anatomically modern humans interbred!
&#10;------------------------------------------------------------------------
&#10;## Bonus content
&#10;As a last sanity check, if we use this time to compute the rate of exponential decay $\lambda$, we should get a nice fit of the theoretical exponential decay curve over the empirical counts of tract lengths in each bin. As a reminder, this is the decay:
&#10;
``` r
plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion of tracts")
```
&#10;Let's try if we can plot the theoretical exponential decay from the estimated time of admixture.
&#10;First, because the exponential tract decay you plotted above shows the decay of the length of introgressed tracts in bins, the $\lambda$ parameter determines the decay expected for tracts falling into to successive bins (so decay of longer and longer tracts). However, our recombination rate $r$ which features in the equation to compute $\lambda$ above describes recombination in units of base pairs, not bins. In order to be able to plot the theoretical exponential decay curve, we therefore first need to compute the average increase in tract lengths from bin to bin:
&#10;
``` r
# compute the average length in each bin
average_bins <- aggregate(length ~ bin, data = tracts, FUN = mean)
&#10;# compute the difference between bins to get "average bin step size"
bin_step <- mean(diff(average_bins$length))
&#10;bin_step
```
&#10;With this, we can overlay the theoretical exponential decay expectation across the empirical distribution of tract lengths in different bins:
&#10;
``` r
r <- 1e-8 # crossovers per bp per generation
t <- 1800 # time of admixture (in generations) we computed above
&#10;# plot the empirical proportions of tracts in each bin
plot(bins, props, xlab = "Neanderthal tract length bin", ylab = "proportion of tracts")
&#10;# overlay theoretical exponential density curve assuming recombination rate r
# and time of introgression t generations ago
lambda <- r * bin_step * t
y <- dexp(bins, rate = lambda)
lines(bins, y, col = "red", lty = 2)
```
&#10;------------------------------------------------------------------------
&#10;If you want to take a closer look at how the tracts data was prepared (it was simulated!), you can see the complete code [here](generate_tracts.R).
-->
