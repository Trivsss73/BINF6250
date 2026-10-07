# Introduction
Gibb's Sampling is a combined Monte Carlo and Markov Chain approach that is useful for motif finding. The approach for motif finding is broken down into three components: initiation, iteration and convergence. We implemented the three components into one function ```GibbsMotifFinder()``` in conjunction with the Driver Program code given to us to produce a plot for the final motif Position Frequency Matrix (PFM). Following the skeleton code, we attempted to construct a code for the ChIP-seq data provided. 

# Pseudocode
Put pseudocode in this box:

```
Notes:
*length of k cannot be longer than the sequence itself
*Pseudocount is established in def build_pwm()
*score the kmer against both forward and reverse complement, and store both separately

Initialization:
For each sequence, pick a random start, extract the kmer, and store as initial motif array

Loop body:
Pick random index i (this is holding out)
build PWM using N-1 motifs (exclude i)
    we get these from previous motifs
use this PWM to go back through the sequences for a new motif using forward and reverse complement
Getting the probabilities:
    currently in log from build pwm, convert back to probabilities
        (* log can be negative, cant have negative probabilities
         * We need probabilities because we are trying to normalize the scores to get the distribution)
    normalizing all the forward and reverse complement scores for each valid position in sequences
replace motifs object with this new sampling for each iteration, maintain strand identity

Convergence:
Use allclose to check convergence for scores from PWM, checking every iteration with a cap of 10,000 if no convergence, storing this round to compare to next round
    #From prof: "until IC doesn't change"

(end loop)
Build final PFM from array


'''
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Kailey: As said in lecture, this definitely was a difficult project- both computationally and conceptually. There were many aspects to consider throughout this project. Gibb's sampling is a new concept to me, so I took the first week to understand the workflow. The lecture slides served as the basis of the pseudocode, with the supplemental lecture aiding in providing additional detail. In addition, the assignment itself had a lot of details so it took time to understand which functions were given and which components that we were required to code ourselves. This project is and was a humbling experience.

# Generative AI Appendix
As per the syllabus
