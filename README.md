## A description of the project

In this project, we implemented the **Smith-Waterman algorithm** for **local pairwise sequence alignment** in Python. 
Pairwise alignment compares two biological sequences (DNA, RNA, or protein) to find regions that are similar, which can point to a functional, structural, or evolutionary relationship. A *local* alignment finds the best-matching **piece** of each sequence, instead of forcing the two whole sequences to line up end to end.

Our program does this in three steps:
1. **Score each cell:** make a scoring matrix with one row per letter of seq1 and one column per letter of seq2, plus an extra first row and column of zeros. Then fill it one cell at a time, row by row, left to right. Each cell looks only at its three neighbors (above-left, above, and left) and takes the best of: the above-left score plus a match or mismatch score, the above score plus a gap penalty, the left score plus a gap penalty, or 0. The 0 is what makes the alignment *local*: a bad stretch resets to 0 instead of dragging the score down. For each cell we also record which option won (the **move**: END, DIAG, UP, or LEFT) in a traceback matrix. If options tie, we pick diagonal first, then up, then left, which prefers pairing letters over adding gaps.

2. **Find the best score:**: after every cell is filled, find the highest score in the matrix. If several cells tie for the highest score, the program picks one of them **at random**; a `seed` setting makes the pick repeatable.

3. **Trace back**: start at the chosen cell and follow the recorded moves backward until reaching a cell whose move is END (a score of 0). Each move adds one column to the front of the alignment.

Match score, mismatch penalty, and gap penalty can all be changed.\
The defaults are match = +1, mismatch = −1, gap = −1.\
We tested the program on the two example pairs of DNA sequences provided in the notebook. We had also added RNA and protein examples. For 1 of the example, we added alignment where scores are different to show that our program allows for customizable gap penalties and match/mismatch scores.

## Pertinent documentation and usage write-up
### Requirements
### Usage
### Project structure

```
project05/
├── project05.ipynb
├── project05.md
├── README.md
```

### Settings
These are the parameters of `smith_waterman()`:
Setting- Default (Meaning)
`seq1`- none (first sequence; its letters label the matrix rows)
`seq2`- none (second sequence; its letters label the matrix columns)
`match`- `1` (added on the diagonal when the two letters are the same)
`mismatch`- `-1` (added on the diagonal when the two letters are different)
`gap`- `-1` (added for an up or left move, when one sequence gets a gap `-`)
`seed`- `None` (random tie-break every run; give a number, e.g. `seed=1`, to repeat the same pick)

### Output
`smith_waterman()` prints:
- `Alignment score: ...` (the highest value in the score matrix)
- `Tied top boxes: N - randomly picked (i, j)`, only when several cells tie for the highest score and returns three things:
- `aligned_seq1`: the aligned part of sequence 1, with `-` where it has a gap
- `aligned_seq2`: the aligned part of sequence 2, with `-` where it has a gap
- `score_matrix`: the filled scoring matrix, of size (len(seq1)+1) × (len(seq2)+1)

### Examples and expected output 
Below we list examples used for testing, note whether or not our output is as expected and and if expected output was provided - note whether or not ours match them.\

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
- For this test, we added alignment where match/mismatch score equal 2 and gap equal = -2
=> result here is also as expected: aligned sequences, different score.
**Test 2: DNA**
- Input provided: seq3 = `GGGGAAAAAAAAAAAACCCC` vs seq4 = `GGGGAAAAAACCCC` 
- Results: see notebook, expected. Two different cells tie for the top score of 10. Both are equally valid local alignments. Our run randomly picked (10, 10). Running with a fixed seed is would make the pick repeatable.

**Test 3: RNA**
- Input (made-up test sequences):seq5 = `AUGGCUAGCA` vs seq6 = `GCUAGCAUGG`
- Results: see notebook, expected. No ties, only one cell has the top score, so this result is the same on every run. 

**Test 4: Protein example**
- Input (Durbin et al., Biological Sequence Analysis) :seq7 = `HEAGAWGHEE` vs seq8 = `PAWHEAE`
- Results: see notebook, three different local alignments tie at the top score of 3. The assumptions is all are equally valid local alignments. Our run randomly picked (10, 10). Running with a fixed seed is would make the pick repeatable.

**Note on the protein test:**
In our program, mismatches cost the same between letter, of which is not the case for protein where amino acids some are very similar to each other (similar size or chemistry), while others are very different. Real tools use substitution matrices like BLOSUM (Blocks Substitution Matrix) or PAM (Point Accepted Mutation) to score how likely one amino acid is to be replaced by another.\
seq7 = `HEAGAWGHEE` vs seq8 = `PAWHEAE` are short made-up-style sequences from a textbook. Real proteins are often hundreds of letters long. So our result shows that the code runs correctly, not necessary about whether or not this is a biologically meaningful alignment. 

**Note on overall limitation of our program:**
- Every gap also costs the same, so for example one long gap costs made out of 3 dashes has the same cost as several short 3 gaps, also 3 dashes; real tools usually charge more to start a gap than to extend it (affine gap penalty), since one long gap is more likely to come from a single event.\
- Our simple scoring creates many ties, and our code shows only one of them.
- Unrelated sequences can match a few letters by chance, and our code doesn’t check whether a score is higher than chance would give.


## Pseudocode
1. Score and direction per cell
```
def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap): #top left -> bottom right  
    #1. Look up the 2 letters (i-1 and j-1 because of the zero border)
    if seq1[i-1] == seq2[j-1]:    
        letter_score = match.      
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
    traceback_matrix = zeros of the same size                         # arrow grid (0 = END)

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

## Reflection narrative that describes what went right, what went wrong, what you didn't understand, how it could have been better, and how you could have been better.
### What went right
### What went wrong
### Things I still don't understand

## Appendix (AI usage, citations, etc.)
**Generative AI version, prompts and justification:**\
Version: Claude (Anthropic),Claude Opus 5.5 Medium\
Prompt: "Give me a make-up RNA seq and protein seq for pairwise alignment algo testing. If there are citation, give it"\
Justification: We wanted to test our algorithm on RNA and protein seq as well. 

**Citations**
Durbin, R., Eddy, S. R., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press.
