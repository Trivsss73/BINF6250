# Introduction
This project introduces Gibbs sampling, a Markov Chain Monte Carlo (MCMC) method
used to discover a shared motif across a set of DNA sequences. The idea is to start
from random motif positions, hold out one sequence at a time, build a position weight
matrix (PWM) from the rest, score every candidate k-mer in the held-out sequence,
and sample a new position probabilistically from the resulting score distribution.
After enough rounds of this, the positions settle on something that actually looks
like a motif, even though no single step is doing any direct optimization.

Sampling is really the whole trick. If we always picked the best-scoring position
each round, the chain would latch onto whatever weak pattern showed up first by
chance and never move. By sampling instead, the low-scoring options still have some
chance of getting picked, which is what lets the chain wander out of a bad spot
and keep looking for a better one.

We ran the sampler in two settings. First on 50 bp promoter regions from Bacillus
subtilis that were pre-filtered to contain AGGAGG (the Shine-Dalgarno motif). Since
we already know the answer here, this one is basically a sanity check that our code
is doing something reasonable. Second on NRF1 ChIP-seq peaks, which is a much bigger
and messier dataset (90,061 sequences of 75 bp, no pre-filter). The second one is
closer to how this algorithm actually gets used in practice, and we are still working
on getting it to converge properly. More on that in the Struggles section.

# Pseudocode
The whole algorithm sits in one function, GibbsMotifFinder. We kept it as a single
block that follows the professor's lecture pseudocode top to bottom, rather than
splitting it into helpers. The three main pieces are initialization, the iteration
loop, and the convergence check.
```
GibbsMotifFinder (seqs, k, seed):
    Inputs:
      seqs - list of DNA strings
      k - motif width
      seed - random seed for reproducibility (default None)
    Output:
      final PFM (4 x k numpy array)

    1. Seed the random number generators.
    2. Guard: if any sequence is shorter than k, raise an error. Also uppercase
       every sequence so case does not break the scoring step.

    Initialization:
    3. For each sequence:
        a. Pick a random start position that still leaves room for a k-mer.
        b. Append that k-mer to the motifs list.
    4. Set up convergence tracking (streak = 0, information_content = None,
       required_streak = 100 as noted in the lecture).

    Iteration loop (up to 10,000 times):
    5. For j in range(10000):
        a. Pick a random sequence index i to be the holdout for this round.
        b. Build a list of every motif except motifs[i].
        c. Build a PFM then a PWM from that list.
        d. Pull out every forward k-mer at every valid start position in seqs[i].
        e. Build the matching reverse complement for each forward k-mer.
        f. Concatenate forward and reverse k-mers into one combined list
           (forward_kmers + reverse_kmers), so no scores get collapsed with max.
        g. Score every k-mer in the combined list against the PWM.
        h. Convert the log2-odds scores to probabilities by exponentiating with
           base 2 (2 ** score), then normalize so they sum to 1.
        i. Sample one index from the probability distribution using a weighted draw.
        j. Replace motifs[i] with the k-mer at the sampled index (already
           reverse-complemented if it came from the reverse half).
        k. Rebuild the full PFM and compute the current IC.
        l. Print the IC every 500 iterations so we can watch the chain climb.

    Convergence check:
    6. If this is not the first iteration and the current IC is within 0.001 bits
       of the previous iteration's IC (np.isclose with atol=1e-3), bump the streak
       up by one. Otherwise reset the streak to zero.
    7. If the streak has reached required_streak (100), break out of the loop.
    8. Save the current IC as the previous_ic for the next round.

    9. Return the final PFM built from the final motifs list.

```

# Successes
Writing the pseudocode first and comparing versions between members made the overall
structure much clearer before anyone wrote actual code. All three of us ended up with
nearly identical algorithms through independent routes, which gave us some confidence
that we were reading the method the same way.

The modular split paid off during debugging. Each helper function could be tested on
its own with a tiny input before being plugged into the main loop. For example, we
tested slide_and_score on a 26 bp fake sequence with AGGAGG planted in the middle,
and watching the score for the correct position dominate the output gave us a lot of
confidence that the scoring step was working before we touched the main algorithm.

A few of the design choices came out of group discussion rather than any one person
figuring them out alone. One was switching from np.exp(score) to 2 ** score for the
probability conversion, since the PWM is already in log2 space and using the matching
base keeps the math honest. Another was keeping the forward and reverse strand scores
separate instead of collapsing them with max. The professor pointed out that
collapsing throws away half the data and amplifies the bias from the AGGAGG pre-filter,
and we built the data structures to keep both strands around the whole way through.

On the Shine-Dalgarno dataset the sampler converged on AGGAGG with final IC around
12 bits after roughly 7000 iterations, matching what we expected to see. The sequence
logo at the end clearly showed the AGG-AGG pattern with weaker flanking positions,
which was a satisfying confirmation that the whole thing worked.

# Struggles
Our first group meeting was rough. All three of us came in with questions and no
one really had a confident picture of what Gibbs sampling was doing or where to
even start the pseudocode. We spent most of that meeting pulling the algorithm
apart out loud, going back and forth on what "probabilistic sampling" meant, why
we were holding out one sequence at a time, and what the PWM was actually for.
By the end we were closer, but we did not have working pseudocode yet. It took
individual reading and a second meeting to get everyone to the same place and
start writing something concrete.

The biggest coding struggle was the convergence check. Our first attempt was
streak-based, meaning we counted how many iterations in a row the IC changed by
less than some tolerance. We started with tolerance 0.01 bits and required streak
of 50. On the Shine-Dalgarno dataset with 837 sequences this triggered "convergence"
after only 51 iterations, with a final IC around 0.47 bits (basically random).
After some digging, we realized the problem. We only update one motif per iteration
out of 837, so the PFM shifts by at most 1/837 per cell per iteration. That means
the per-iteration IC change is naturally way below 0.01 bits whether or not the
chain is actually settling. The streak was triggering on the natural slowness of
the chain, not on real stability. We tightened the tolerance to 0.001 bits, which
fixed the Shine-Dalgarno case and let the chain run long enough to actually find
AGGAGG.

NRF1 brought the same problem back on a bigger scale. With 90,061 sequences, the
per-iteration change is naturally below 0.001 bits too, and the chain was still
converging prematurely with final IC near 0.68 bits. We added a min_iter floor to
the convergence check so the chain is forced to run for a minimum number of
iterations before any early stopping is allowed. Even with this, picking the right
combination of tol, required_streak, and min_iter for NRF1 is still an open question
for us. The chain is slow, the signal is weaker without a pre-filter, and we are
still running longer experiments to see if it converges on anything meaningful.

The professor was clear in lecture that we should not use max() to collapse the
forward and reverse strand scores at each position, since that cuts the data in
half and amplifies the AGGAGG pre-filter bias. We kept the scores separate the
whole way through to avoid this.

Some smaller stumbling blocks came up along the way. The off-by-one math for valid
start positions (len(seq) - k + 1 is the number of valid starts, not an upper bound)
tripped us up at first. The interleaved score list meant we had to be careful
recovering the position (idx // 2) and strand (strand_labels[idx]) from a single
sampled index. And it was easy to forget that the holdout sequence being scored is
the full-length seqs[i], not the k-mer motifs[i] we had been tracking in the main
state.

Towards the end we struggled with Github as a team, once we accepted the changes, the pull request automatically closed, we would love inputs on how to fix or deal with this.

# Personal Reflections
## Group Leader
Dhaivat - Leading this project for the first time was a big learning curve. Early
on, when all three of us came into the first meeting with questions and no real
grip on the algorithm, it felt like the answer to every question we asked was just
"it's random." That was frustrating because it did not give us anything concrete
to code against. It took watching the professor's walkthrough video and sitting
with the lecture flowchart to finally see the structure under all the randomness,
although I will admit I misread parts of the flowchart on my first pass. The
"holdout PWM" label specifically threw me because I kept thinking it meant a PWM
built from the holdout sequence, when it actually means the PWM we use to score
the holdout. That mix-up cost me some time before I realized what was going on.

The algorithm itself took time to click. What finally helped was realizing we
are sampling, not optimizing, and that the randomness is actually the point.
Some parts still feel fuzzy to me though. The deeper theory of why Gibbs sampling
converges to the right distribution is something I can follow at the surface but
not really explain, and I would want to come back to that at some point.

The convergence check was the most educational struggle. My original plan was a
window-based check that looked at the range of IC values over the last N iterations,
but Maranda had already written a streak-based version that matched the professor's
"until IC doesn't change" phrasing more closely, so we went with hers. That turned
out to have a subtle bug on big datasets. Our first run told us the sampler had
converged after 51 iterations with IC near 0.47 bits, which was nowhere near a
real motif. Working through why (one update per iteration out of 837 sequences
means the per-iteration IC change is naturally tiny) and adding a min_iter floor
was more satisfying than just tweaking the tolerance would have been.

NRF1 has been genuinely annoying. The sampler runs for thousands of iterations
and the IC barely moves off the noise floor, so the sequence logo plot ends up
blank. We know the dataset is bigger and noisier than Shine-Dalgarno and that
convergence is going to be slower, but after multiple parameter tweaks it still
has not produced anything that looks like a motif. I do not have a good answer
for why yet, and that is still bothering me.

Setting up the group repo this time was also new. Managing branches, PRs, and
pulling in teammate work from the owner side was different from being a contributor,
and I leaned on my group members to figure out some of the git workflow pieces.

## Other member
Kailey: As said in lecture, this definitely was a difficult project- both computationally and conceptually. There were many aspects to consider throughout this project. Gibb's sampling is a new concept to me, so I took the first week to understand the workflow. The lecture slides served as the basis of the pseudocode, with the supplemental lecture aiding in providing additional detail. In addition, the assignment itself had a lot of details so it took time to understand which functions were given and which components that we were required to code ourselves. This project is and was a humbling experience. 

The logic of the project stemmed from me creating a concept map of Gibb's sampling. From there I was able to construct the pseudocode as broken down into the 3 main steps: initialization, iteration and convergence. Implementing the actual code, is the most difficult to me. Our group met a total of three times to discuss the project- first to overview the material, second to assess pseudocode and overview of the supplemental lectures and third to review our progress. Our group selected an independent approach when it came to code and would collaborate to view which aspects of each others code we liked to create one cohesive code. The ChIP-seq challenge, had my brain turning as to figure out what could the output be. From this project, I learned that it is okay not to complete the coding, but to put our best effort to understand the concept and reasoning behind this specific algorithm. 

Mara: This was quite a difficult project, with many new concepts to understand before implementation, even with the extra week to complete it. Working with my group mates was very beneficial to brainstorm how the project should be implemented, and even still there was a lot of debugging. Overall I think this project was very useful to start understanding some more difficult concepts that I will eventually encounter in this field, improving my collaboration with others, and still getting more comfortable with GitHub and Jupyter Notebooks. 

# Generative AI Appendix
Claude Opus 4.7 was used to understand the algorithm, break it down in parts and undestand how it
applies to MOTIF finding. It was also used to understand the possible ways to deal with 
convergence. Lastly to create debugging code to make sure the functions are working as designed. 

Prompt: Given the pseudocode, provide step by step how I can approach the technical coding.
Output: Claude provided explanation and reasoning on how to construct the code without giving the answers.
Reasoning: With the complexity of this project, I wanted to create a clear step-by-step process on how to improve my technical coding experience and have an organized to-do list.
 
Prompt: Given the code, provide feedback on how to debug the following errors.
Output: Claude provided reasonings to my errors, such as Value Errors and Import Errors (because I forgot to run the import block first).
Reasoning: In one of my courses, it was mentioned that after about 20 minutes of struggling with debugging, it is okay to ask generative AI to assist. Here, I was struggling a bit, just for a few simple fixes to my code.  

Prompt: I am working on a Gibb's Sampling project, I am required to develop the algorithm from scratch but first I would like to understand it first. How this algorithm is used in bioinformatics. explain the algorithm clearly, not using too much technical jargon.
Output: Claude explained the algorithm, the stats behind it, how it applies to bioinformatics. 

Prompt: Explained the output of the project, gave a draft pseudocode and asked for improvements, point out the gaps in understanding and workflow. 
Output: gave its opinions on what works and why?, what would not work? and gave possible ways to tryout.
