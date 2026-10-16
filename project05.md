Project 05: Pairwise Alignments
================
20251008

- [Learning Objectives](#learning-objectives)
- [Background](#background)
- [Smith-Waterman Algorithm Implementation
  Activity](#smith-waterman-algorithm-implementation-activity)
  - [Introduction](#introduction)
  - [Activity Overview](#activity-overview)
  - [Instructions](#instructions)
  - [Follow these steps:](#follow-these-steps)
  - [Imports](#imports)
  - [Implement Smith-Waterman
    algorithm](#implement-smith-waterman-algorithm)

## Learning Objectives

1.  Conceptually understand dynamic programming and sequence alignment
2.  Implement Smith-Waterman algorithm for local alignment

------------------------------------------------------------------------

## Background

This week we will be implementing Smith-Waterman. This is a dynamic
programming algorithm used for local sequence alignment.

As a reminder, the scoring for Smith-Waterman only uses the scores from
the positions above, left, and above-left of the current position in the
matrix as below:

<center>
<img src="figures/Smith-Waterman_scoring.png">
</center>

For traceback, you will need to keep track of the direction of the
arrows in a matrix and then begin traceback from the maximum value.

# Smith-Waterman Algorithm Implementation Activity

## Introduction

Pairwise sequence alignment is a fundamental technique in bioinformatics
used to compare two biological sequences, such as DNA, RNA, or proteins.
It helps identify regions of similarity that may indicate functional,
structural, or evolutionary relationships between the sequences.

The brute force method for sequence alignment, which involves comparing
every possible alignment, is highly inefficient. For two sequences of
length m and n, the time complexity would be
$T(m,n) = O(mn \cdot 2^{m+n})$, making it impractical for longer
sequences.

Dynamic programming offers a more efficient solution. It breaks down the
problem into smaller subproblems and stores their solutions to avoid
redundant computations. This approach reduces the time complexity to
$O(mn)$.

The Smith-Waterman algorithm is a dynamic programming approach for local
sequence alignment. It identifies the optimal local alignment between
two sequences by comparing all possible pairs of segments from the
sequences and finding the best-scoring alignment.

## Activity Overview

In this activity, you will implement the Smith-Waterman algorithm for
local sequence alignment. The algorithm consists of two main parts:

Scoring matrix calculation Traceback for optimal alignment
reconstruction

The pseudocode for the Smith-Waterman algorithm is as follows:

$$
\begin{aligned}
& \text{Initialize scoring matrix } H \text{ with zeros} \\
& \text{For } i \text{ from } 1 \text{ to } m: \\
&     \text{For } j \text{ from } 1 \text{ to } n: \\
&         H[i][j] = \max \begin{cases}
&             0 \\
&             H[i-1][j-1] + \text{match/mismatch score} \\
&             H[i-1][j] + \text{gap penalty} \\
&             H[i][j-1] + \text{gap penalty}
&         \end{cases} \\
& \text{Find the highest score in } H \text{ and its position} \\
& \text{Perform traceback from the highest score position} \\
& \text{Return the optimal local alignment}
\end{aligned}
$$

## Instructions

Working in pairs, implement the Smith-Waterman algorithm using Python.
Your implementation should include:

1.  A function to calculate the scoring matrix (`cal_score`)
2.  A function to perform the traceback and reconstruct the optimal
    alignment (`traceback`)
3.  A dynamic scoring strategy that allows for customizable gap
    penalties and match/mismatch scores (`smith_waterman`)

## Follow these steps:

### Implement the scoring matrix calculation function:

1.  Initialize the matrix with zeros
2.  Fill the matrix using the Smith-Waterman algorithm
3.  Return the completed matrix and the position of the highest score

### Implement the traceback function

1.  Start from the highest score position
2.  Trace back through the matrix to reconstruct the optimal alignment
3.  Return the aligned sequences and the alignment score

### Implement a main function that:

1.  Takes two sequences as input
2.  Accepts parameters for match score, mismatch penalty, and gap
    penalty
3.  Calls the scoring and traceback functions
4.  Prints the aligned sequences and the alignment score

Test your implementation with various input sequences and scoring
parameters

------------------------------------------------------------------------

## Imports

``` python
import numpy as np
```

------------------------------------------------------------------------

## Implement Smith-Waterman algorithm

    SmithWaterman(seq1, seq2, match, mismatch, gap)
        Initialize [len(seq1)+1] x [len(seq2)+1] numpy array as scoring matrix with first column and row equal to 0
        
        Fill scoring matrix (score_matrix) and Traceback matrix, record the position with max score (max_pos):
        for i in each row number:
            for j in each column number:
                S[i][j] = max( S[i-1][j-1] + compute_diag_score, S[i-1][j] + gap_score, H[i][j-1] + gap_score, 0 )
                T[i][j] = direction of max( S[i-1][j-1] + compute_diag_score, S[i-1][j] + gap_score, S[i][j-1] + gap_score, 0 )
        
        Traceback. Find the optimal path through scoring matrix starting at max_pos
        

``` python
def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):
    '''Calculate score for position (i,j) in scoring matrix, also record move to trace back
    
    Args:
        matrix (numpy array): scoring matrix
        seq1 (str): sequence 1
        seq2 (str): sequence 2
        i (int): current row number
        j (int): current column number
        
    Returns:
        score in position (i,j)    
        move to trace back: 0-END, 1-DIAG, 2-UP, 3-LEFT
        
    Pseudocode:
        Calculate scores based on upper-left, up, and left neighbors:
            diag_score = upper-left + (match or mismatch)
            up_score = up + gap
            left_score = left + gap
        score = max(0, diag_score, up_score, left_score)
        traceback = maximum direction or end
        
    '''
    pass
    
```

``` python
def traceback(seq1, seq2, traceback_matrix, maximum_position):
    '''Find the optimal path through scoring marix
        
        diagonal: match/mismatch
        up: gap in seq1
        left: gap in seq2
        
    Args:
        seq1 (str) : First sequence being aligned
        seq2 (str) : Second sequence being aligned
        traceback_matrix (numpy array): traceback matrix
        maximum_position (tuple): starting position to trace back from
        
    Returns:
        aligned_seq1 (str): e.g. GTTGAC
        aligned_seq2 (str): e.g. GTT-AC
        
    Pseudocode:
        while current_move != END:
            current_move = traceback_matrix[current_row][current_col]
            if current_move == DIAG:
                ...
            elif current_move == UP:
                ...
            elif current_move == LEFT:
                ...
            
    '''
    pass
            
```

``` python
def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
    '''Smith-Waterman algorithm for local alignment
    
    Args:
        seq1 (str): input seq 1
        seq2 (str): input seq 2
        match: default = +1
        mismatch: default = -1
        gap: default = -1
    
    Returns:
        aligned_seq1 (str)
        aligned_seq2 (str)
        score_matrix (numpy array): scoring matrix
    '''
    pass
```
