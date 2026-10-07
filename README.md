# Introduction
Gibb's Sampling is a combined Monte Carlo and Markov Chain approach that is useful for motif finding. The approach for motif finding is broken down into three components: initiation, iteration and convergence. We implemented the three components into one function ```GibbsMotifFinder()``` in conjunction with the Driver Program code given to us to produce a plot for the final motif Position Frequency Matrix (PFM). Following the skeleton code, we attempted to construct a code for the ChIP-seq data provided. 
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
The main driver is one function (GibbsMotifFinder) that runs the convergence loop.
We broke the smaller operations into their own helper functions so we could test
each piece on its own before plugging them into the loop.
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
GibbsMotifFinder(seqs, k, seed, max_iter, required_streak, tol, min_iter):
    Inputs:
      seqs - list of DNA strings
      k - motif width
      seed - random seed for reproducibility
      max_iter - hard cap on iterations
      required_streak - how many stable iterations before calling it converged
      tol - how small the IC change has to be to count as stable
      min_iter - floor on iterations before convergence can trigger
    Output:
      final PFM (4 x k numpy array)

    1. Seed the random generators.
    2. Guard: if any sequence is shorter than k, bail out with an error.
    3. Initialize motifs, strands, and positions by calling init_motifs.
    4. Set up the convergence tracking state (streak = 0, previous_ic = None).
    5. Loop up to max_iter times:
        a. Pick one random sequence index i to be the holdout for this round.
        b. Build a list of every motif except motifs[i].
        c. Build a PWM from that list using init_pwm.
        d. Slide across the full-length seqs[i] and score every candidate k-mer
           on both strands using slide_and_score.
        e. Sample a new motif, strand, and position from the score distribution
           using sample_new_motif.
        f. Replace motifs[i], strands[i], positions[i] with the new picks.
        g. Rebuild the full PFM, compute IC, and update the convergence streak.
        h. If the streak hit required_streak AND we are past min_iter, break.
    6. Return the final PFM built from the current motifs.


init_motifs (seqs, k, rng):
    Inputs:
      seqs - list of DNA strings
      k - motif width
      rng - numpy random generator
    Output:
      motifs, strands, positions - three parallel lists, one entry per sequence

    1. For each sequence:
        a. Pick a random start position that still leaves room for a k-mer.
        b. Extract the k-mer at that position.
        c. Flip a coin for the strand. If reverse, replace the k-mer with its
           reverse complement and label it 'R'. Otherwise label it 'F'.
    2. Return the three lists.


init_pwm (other_motifs, k):
    Inputs:
      other_motifs - the k-mers from every sequence except the holdout
      k - motif width
    Output:
      pwm - 4 x k log2 odds matrix

    1. Build a PFM from other_motifs.
    2. Convert it to a PWM.
    3. Return the PWM.


slide_and_score (holdout_seq, pwm, k):
    Inputs:
      holdout_seq - the full-length sequence we are scoring
      pwm - the PWM to score against
      k - motif width
    Output:
      scores and strand_labels - two parallel lists, interleaved F and R per position

    1. For each valid start position in the holdout sequence:
        a. Pull out the forward k-mer and its reverse complement.
        b. Score both against the PWM.
        c. Append the forward score and 'F', then the reverse score and 'R',
           so position p lands at indices 2p and 2p+1.
    2. Return the two lists.


sample_new_motif (scores, strand_labels, holdout_seq, k, rng):
    Inputs:
      scores - the full interleaved F/R score list
      strand_labels - matching strand labels
      holdout_seq - the sequence we are sampling from
      k - motif width
      rng - numpy random generator
    Output:
      the new k-mer, its strand, and its position in the holdout

    1. Convert the scores to probabilities. Since the PWM is in log2 space, we
       exponentiate with base 2 (2 ** score), then normalize so they sum to 1.
    2. Pick one index from the probability distribution using a weighted draw.
    3. Recover the position (idx // 2) and the strand (strand_labels at idx).
    4. Extract the k-mer at that position from the holdout sequence.
    5. If the strand is 'R', reverse complement the k-mer.
    6. Return the k-mer, strand, and position.


check_convergence (current_ic, previous_ic, streak, required_streak, tol):
    Inputs:
      current_ic - IC after this iteration
      previous_ic - IC from the last iteration, or None on the first call
      streak - current streak count going in
      required_streak - how long the streak has to be before we call it converged
      tol - how small the IC change has to be to count as stable
    Output:
      updated streak and a converged flag (True / False)

    1. If this is the first iteration (previous_ic is None), reset streak to 0
       and return False.
    2. If the IC change is smaller than tol, bump the streak up by one.
    3. Otherwise, reset the streak to zero.
    4. Return the streak and whether it has hit required_streak.

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
