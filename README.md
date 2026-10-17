## A description of the project

In this project, we implemented the **Smith-Waterman algorithm** for **local pairwise sequence alignment** in Python. 
Pairwise alignment compares two biological sequences (DNA, RNA, or protein) to find regions that are similar, which can point to a functional, structural, or evolutionary relationship. A *local* alignment finds the best-matching **piece** of each sequence, instead of forcing the two whole sequences to line up end to end.

Our program does this in three steps:
1. **Score each cell:** fill a scoring matrix one cell at a time. Each cell looks only at its three neighbors (above-left, above, and left) and takes the best of: the above-left score plus a match or mismatch score, the above score plus a gap penalty, the left score plus a gap penalty, or 0. The 0 is what makes the alignment *local*: a bad stretch resets to 0 instead of dragging the score down. For each cell we also record which direction the best score came from (the **move**) in a traceback matrix.
2. **Find the best score:** keep track of the highest score in the matrix and where it is.
3. **Trace back:** start at the highest score and follow the recorded moves backward until reachingg a cell whose move is END (a score of 0). Each move adds one column to the alignment.

Match score, mismatch penalty, and gap penalty can all be changed.\
The defaults are match = +1, mismatch = −1, gap = −1.\
We tested the program on the two example pairs of DNA sequences provided in the notebook. There is no external data set.
=>>>>. add 1 more example for RNA and 1 example for protein.

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
### Output

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
def smith_waterman(seq1, seq2, match = 1, mismatch = -1, gap = -1):
    #1. Set up 2 empty grids
    score_matrix     = zeros of size (len(seq1)+1) x (len(seq2)+1) #number grid
    traceback_matrix = zeros of the same size #direction grid
    #2. For each cell, row by row, left to right, update the score:
    # Set up sticky note for currently the best number and direction (0):
    best_score = 0, best_position = (0, 0) 
    # Call on function cal_score to hand back the number and arrow to the cell, update the note with it
    for i from 1 to len(seq1):
        for j from 1 to len(seq2):
            score_matrix[i][j], traceback_matrix[i][j] = cal_score(score_matrix, seq1, seq2, i, j, match, mismatch, gap)
            if score_matrix[i][j] > best_score:
=>>>>>>>>>> Marcus said using random instead of 1st pick !!!!!!!

                best_score = score_matrix[i][j]
                best_position = (i, j)
    # Call on traceback function to get aligned seq 1 and 2
    aligned_seq1, aligned_seq2 = traceback(seq1, seq2, traceback_matrix, best_position)
    return aligned_seq1, aligned_seq2, score_matrix
```
## Reflection narrative that describes what went right, what went wrong, what you didn't understand, how it could have been better, and how you could have been better.
### What went right
### What went wrong
### Things I still don't understand

## Appendix (AI usage, citations, etc.)
def traceback(seq1, seq2, traceback_matrix, maximum_position):
