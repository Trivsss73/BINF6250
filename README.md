# Introduction
In this project, we implemented the **Smith-Waterman algorithm** for **local pairwise sequence alignment** in Python. 
Pairwise alignment compares two biological sequences (DNA, RNA, or protein) to find regions that are similar, which can point to a functional, structural, or evolutionary relationship. A *local* alignment finds the best-matching **piece** of each sequence, instead of forcing the two whole sequences to line up end to end.

Our program does this in three steps:
1. **Score each cell:** make a scoring matrix with one row per letter of seq1 and one column per letter of seq2, plus an extra first row and column of zeros. Then fill it one cell at a time, row by row, left to right. Each cell looks only at its three neighbors (above-left, above, and left) and takes the best of: the above-left score plus a match or mismatch score, the above score plus a gap penalty, the left score plus a gap penalty, or 0. The 0 is what makes the alignment *local*: a bad stretch resets to 0 instead of dragging the score down. For each cell we also record which option won (the **move**: END, DIAG, UP, or LEFT) in a traceback matrix. If options tie, we pick diagonal first, then up, then left, which prefers pairing letters over adding gaps.

2. **Find the best score:** after every cell is filled, find the highest score in the matrix. If several cells tie for the highest score, the program picks one of them **at random**; a `seed` setting makes the pick repeatable.

3. **Trace back**: start at the chosen cell and follow the recorded moves backward until reaching a cell whose move is END (a score of 0). Each move adds one column to the front of the alignment.

Match score, mismatch penalty, and gap penalty can all be changed.\
The defaults are match = +1, mismatch = −1, gap = −1.\
We tested the program on the two example pairs of DNA sequences provided in the notebook. We had also added RNA and protein examples. For 1 of the example, we added alignment where scores are different to show that our program allows for customizable gap penalties and match/mismatch scores.

# Documentation and Results
### Requirements
- Python 3 (Python 3.12.4)
- Jupyter Notebook, we used PyCharm
- numpy, used for the score and traceback matrices and for the random tie-break
  
### Usage
Open project05.ipynb in Jupyter and run all cells in order, top to bottom. \
The notebook is split into:
Imports: import numpy as np and the move codes END, DIAG, UP, LEFT = 0, 1, 2, 3\
Core functions: cal_score() scores one cell and records its move, traceback() rebuilds the alignment from the traceback matrix, and smith_waterman() is the main function that fills both matrices, picks the top-scoring cell, and calls traceback().\
Tests: Test 1 (DNA, provided, plus a custom-scoring run), Test 2 (DNA, provided), Test 3 (RNA), and Test 4 (protein)

### Project structure

```
project05/
├── project05.ipynb
├── project05.md
├── README.md
```

### Settings
These are the parameters of `smith_waterman()`:\
Setting- Default (Meaning)\
`seq1`- none (first sequence; its letters label the matrix rows)\
`seq2`- none (second sequence; its letters label the matrix columns)\
`match`- `1` (added on the diagonal when the two letters are the same)\
`mismatch`- `-1` (added on the diagonal when the two letters are different)\
`gap`- `-1` (added for an up or left move, when one sequence gets a gap `-`)\
`seed`- `None` (random tie-break every run; give a number, e.g. `seed=1`, to repeat the same pick)

### Output
`smith_waterman()` prints:
- `Alignment score: ...` (the highest value in the score matrix)
- `Tied top boxes: N - randomly picked (i, j)`, only when several cells tie for the highest score and returns three things:
- `aligned_seq1`: the aligned part of sequence 1, with `-` where it has a gap
- `aligned_seq2`: the aligned part of sequence 2, with `-` where it has a gap
- `score_matrix`: the filled scoring matrix, of size (len(seq1)+1) × (len(seq2)+1)

### Examples and outputs 
Below we list examples used for testing, note whether or not our output is as expected and if expected output was provided - note whether or not ours match them.

**Test 1: DNA**
- Input provided: seq1 = `TACTTAG` vs seq2 = `CACATTAA` 
- Expected output provided (ours matched it):
```
AC-TTA
ACATTA
[[0 0 0 0 0 0 0 0 0]
 [0 0 0 0 0 1 1 0 0]
 [0 0 1 0 1 0 0 2 1]
 [0 1 0 2 1 0 0 1 1]
 [0 0 0 1 1 2 1 0 0]
 [0 0 0 0 0 2 3 2 1]
 [0 0 1 0 1 1 2 4 3]
 [0 0 0 0 0 0 1 3 3]]
```
Only one cell has the top score, so the result is the same every run. 
- For this test, we added alignment where match score equal 2 and gap equal = -2 while letting mismatch stay -1.
=> result here is also as expected: aligned sequences, different score.

**Test 2: DNA**
- Input provided: seq3 = `GGGGAAAAAAAAAAAACCCC` vs seq4 = `GGGGAAAAAACCCC` 
- Results: see notebook, expected. Two different cells tie for the top score of 10. Both are equally valid local alignments. Our run randomly picked (20, 14). Running with a fixed seed would make the pick repeatable.

**Test 3: RNA**
- Input (made-up test sequences):seq5 = `AUGGCUAGCA` vs seq6 = `GCUAGCAUGG`
- Results: see notebook, expected. No ties, only one cell has the top score, so this result is the same on every run. 

**Test 4: Protein example**
- Input (Durbin et al., Biological Sequence Analysis) :seq7 = `HEAGAWGHEE` vs seq8 = `PAWHEAE`
- Results: see notebook, three different local alignments tie at the top score of 3. The assumption is all are equally valid local alignments. Our run randomly picked one of the tied cells, which can change between runs since the seed is not fixed. Running with a fixed seed would make the pick repeatable.

**Note on the protein test:**
In our program, mismatches cost the same between letter, of which is not the case for protein where amino acids some are very similar to each other (similar size or chemistry), while others are very different. Real tools use substitution matrices like BLOSUM (Blocks Substitution Matrix) or PAM (Point Accepted Mutation) to score how likely one amino acid is to be replaced by another.\
seq7 = `HEAGAWGHEE` vs seq8 = `PAWHEAE` are short made-up-style sequences from a textbook. Real proteins are often hundreds of letters long. So our result shows that the code runs correctly, not necessarily about whether or not this is a biologically meaningful alignment.

# Pseudocode
1. Score and direction per cell
```
def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap): #top left -> bottom right  
    #1. Look up the 2 letters (i-1 and j-1 because of the zero border)
    if seq1[i-1] == seq2[j-1]:    
        letter_score = match      
    else:
        letter_score = mismatch

    #2. Calculate four score options (diagonal, down, sideway, 0), choose the biggest 
    diag_score = matrix[i-1][j-1] + letter_score 
    up_score   = matrix[i-1][j]   + gap
    left_score = matrix[i][j-1]   + gap
    score = max(0, diag_score, up_score, left_score)

    #3. For each chosen option (recorded score per cell), note the the direction it was calculated from (if tie, then diagonal first, then down, then sideways)
    if score == 0:               move = END
    else if score == diag_score: move = DIAG
    else if score == up_score:   move = UP
    else:                        move = LEFT
    # Return score and direction
    return score, move
```
2. Optimal path via traceback 
```
def traceback(seq1, seq2, traceback_matrix, maximum_position): # right to left
    aligned_seq1 = "", aligned_seq2 = ""
    # Start at cell with the biggest score, going over the noted directions until landing on "END":
    i, j = maximum_position
    while traceback_matrix[i][j] != END:
        #1. If the direction is up left diagonal => add to seq 1 and seq 2 the aligned / matching row / column letter
        move = traceback_matrix[i][j]
        if move == DIAG:
            add seq1[i-1] to the front of aligned_seq1
            add seq2[j-1] to the front of aligned_seq2
            i = i - 1, j = j - 1

        #2. If not, it's a gap, add "-" 
        else if move == UP:
            add seq1[i-1] to the front of aligned_seq1
            add "-" to the front of aligned_seq2
            i = i - 1
        else if move == LEFT:
            add "-" to the front of aligned_seq1
            add seq2[j-1] to the front of aligned_seq2
            j = j - 1
    return aligned_seq1, aligned_seq2
```
3. Main function: Smith-Waterman algorithm for local alignment
```
def smith_waterman(seq1, seq2, match = 1, mismatch = -1, gap = -1, seed = None):
    # 1. Set up the random picker (only used if several boxes tie for the top score)
    rng = random number maker using seed

    # 2. Set up 2 empty grids
    score_matrix     = zeros of size (len(seq1)+1) x (len(seq2)+1)   # number grid
    traceback_matrix = zeros of the same size                        # direction grid (0 = END)

    # 3. Fill every box, row by row, left to right
    #    cal_score hands back the number and the arrow for each box
    for i from 1 to len(seq1):
        for j from 1 to len(seq2):
            score_matrix[i][j], traceback_matrix[i][j] = cal_score(score_matrix, seq1, seq2, i, j, match, mismatch, gap)

    # 4. After the loops: find the top score and pick where to start the traceback
    best_score     = the biggest number in score_matrix
    tied_positions = every (i, j) where score_matrix equals best_score
    best_position  = pick one of tied_positions at random

    # 5. Call the traceback function to get aligned seq 1 and 2
    aligned_seq1, aligned_seq2 = traceback(seq1, seq2, traceback_matrix, best_position)

    # 6. Print the score (and the tie, if any), then return
    print best_score
    if more than one tied position: print how many tied and which one was picked
    return aligned_seq1, aligned_seq2, score_matrix
```

# Successes
Before either of us wrote any code, we spent a good amount of time discussing how dynamic programming actually works and how it applies to Smith-Waterman. We talked through the idea of breaking a big problem into smaller overlapping subproblems and storing each answer so it never has to be recomputed. For alignment, that means every cell in the matrix is the best score for a pair of prefixes, built only from the three neighbors that are already filled in. Getting on the same page about that first made the scoring and traceback steps feel much more natural once we started coding, instead of just following the pseudocode line by line.

After that, we each wrote our own version of the algorithm, then compared them side by side before picking one for submission. We also each came up with our own way to settle ties when more than one cell shares the top score. Dhaivat tracked the max position during the fill loop with `>=`, which keeps the last tied cell it finds, while Trang collected every tied cell after the fill and picked one at random, with a `seed` to make the pick repeatable. Both versions produced the exact expected matrix and alignment (`AC-TTA` / `ACATTA`) on the slide example, which gave us confidence the scoring and traceback logic were right before we started making design choices.

Comparing the two versions is where most of the useful discussion happened. Besides ties, they differed in a few other small ways (how the direction was picked in `cal_score` and how the aligned strings were built), and talking through each one made us think about why it mattered instead of just whether the code ran. We went with Trang's version because it finds every cell tied for the top score after the fill, reports how many there are, and uses a `seed` so the pick can be repeated.

Hand-tracing the traceback also helped us catch that the template docstring had UP and LEFT swapped. With seq1 on the rows, an UP move puts the gap in seq2 and a LEFT move puts the gap in seq1, and our docstrings now say that.

On top of the two provided tests, we added a custom-scoring run (`match=2, gap=-2`), an RNA pair, and a protein pair, to cover the "test with various inputs and scoring parameters" part of the assignment.

# Struggles
The biggest design question was how to handle ties. On test 2, two different cells tie at a score of 10 (`GGGGAAAAAA` and `AAAAAACCCC`), and both are equally valid local alignments. A version that quietly keeps one of them gives a reproducible answer but hides that a tie existed, while a random pick is more honest but can change between runs. We settled on reporting the tie count and using a `seed`, but we're still not sure it's the best option compared to just returning every tied alignment.

**Note on overall limitation of our program:**
- Every gap also costs the same, so for example one long gap costs made out of 3 dashes has the same cost as several short 3 gaps, also 3 dashes; real tools usually charge more to start a gap than to extend it (affine gap penalty), since one long gap is more likely to come from a single event.\
- Our simple scoring creates many ties, and our code shows only one of them.\
- Unrelated sequences can match a few letters by chance, and our code doesn’t check whether a score is higher than chance would give.

**Things we still don't understand:**
- How substitution matrices like BLOSUM are built and why some amino-acid swaps score higher.
- How affine gap penalties change the scoring step.
- Whether a random pick is the best way to handle ties, or whether reporting all tied alignments would be better, or whether and what is the reason random pick is better than in the list of cells with highest score, pick the earliest encounter in the walk.

# Personal Reflections
## Group Leader
Dhaivat- 

Breaking the problem down before writing any code helped a lot. Before touching `cal_score`, I worked out by hand which two letters each cell compares, and figured out that cell `(i, j)` has to look at `seq1[i-1]` and `seq2[j-1]` because of the extra zero row and column. I prototyped it by scoring a couple of single cells on paper, like `(2, 2)` and `(1, 1)`, and checking them against the expected matrix before writing the loops. When I finally ran the full matrix, it matched the expected output on the first try, which honestly felt great.

I also wrote my own version of the full algorithm before we compared notebooks. A few things in mine were different from Trang's:
- In `cal_score` I used `np.argmax` on the list `[0, diag, up, left]`. Since argmax returns the first max, the order of the list itself becomes the tie-break rule (END, then DIAG, then UP, then LEFT).
- In `traceback` I added letters to the end of the strings and flipped both at the very end with `[::-1]`, instead of adding to the front each step.
- I tracked the max position during the fill loop using `>=`, so if two cells tied for the top score, my code kept the last one it found.

Hand-tracing the traceback also helped me catch that the template docstring had UP and LEFT swapped. With seq1 on the rows, a LEFT move puts the gap in seq1, not seq2.

Honestly, this assignment felt a lot more relaxed than the last two. After de Bruijn and Gibbs sampling, where we spent a lot of time debugging convergence, recursion limits, and huge datasets, Smith-Waterman was much easier to wrap my head around. Most of what went wrong this time was small stuff.

One silly mistake I made and didn't catch at first glance was writing `pair_score = -mismatch`. Since `mismatch` is already passed in as -1, that turned every mismatch into +1. It took me a minute to see that the negative sign lives in the value, not in my code, and the same goes for `gap`.

I also kept mixing up `i` and `j` with the two sequences in traceback (writing `seq1[j-1]` at one point). I had to remind myself that seq1 goes with rows and `i`, and seq2 goes with columns and `j`.

The `>=` choice ended up mattering more than I expected. On test 2, there are two cells tied at 10, and my version always returned `AAAAAACCCC`, while Trang's random pick could return either one. When we compared, I realized my version silently picked one answer without telling the user a tie existed. That's a big reason we went with Trang's approach for the final submission: it finds every tied cell after the fill and reports the tie count, with a `seed` for repeatability.

## Other members
Trang -
Overall, this project feels more relaxed. It was very helpful to talk through and align understanding of the algorithm and pseudo code prior to writing anything. For ties, my first version kept the first top-scoring cell it found (strict >), similar to Trivedi. I changed it to collect every tied cell after the fill and pick one at random, with a seed so a run can be repeated. On the protein test, three alignments tie at 3, so the output can change between runs. We discussed some of the limitation of our program both overall and for protein specifically. A next step for the program would be reporting every tied alignment instead of one, and adding a substitution matrix and affine gap penalties for proteins.

# Generative AI Appendix
**Generative AI:**\
We used Claude (Anthropic) to: explain and check our understanding of the Smith-Waterman algorithm; suggest the RNA and protein test sequences and cite where it was taken from, or not if just made up. 

**Citations**
Durbin, R., Eddy, S. R., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press.
