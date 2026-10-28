"""Week 5 Smith-Waterman local alignment helpers, reused in Project 06.

Direction codes used in the traceback matrix:
    0 END, 1 DIAG, 2 UP, 3 LEFT

Notes:
    - Score matrices are int, so only whole-number scoring values are supported.
      Project 06 always uses match 2, mismatch -1, gap -1 (see README).
    - The default scores below (1 / -1 / -1) are the week 5 defaults. Project 06
      always passes its own values in explicitly, so these defaults are never used there.
"""

import numpy as np


def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):
    """Scores one cell of the Smith-Waterman matrix.

    Args:
        matrix (np.ndarray): score matrix filled up to (i, j)
        seq1 (str): sequence along the rows
        seq2 (str): sequence along the columns
        i (int): row of the cell (1-based, row 0 is padding)
        j (int): column of the cell (1-based, column 0 is padding)
        match (float): score for matching characters
        mismatch (float): penalty for mismatched characters (negative)
        gap (float): penalty for a gap (negative)

    Returns:
        tuple: (best score for this cell, direction code 0-3)
    """
    # i-1 and j-1 because row 0 and col 0 are just padding
    if seq1[i-1] == seq2[j-1]:
        pair_score = match
    else:
        pair_score = mismatch   # mismatch is already negative, no minus sign needed

    diag_score = matrix[i-1][j-1] + pair_score
    up_score   = matrix[i-1][j] + gap
    left_score = matrix[i][j-1] + gap

    # list position = direction code (0 END, 1 DIAG, 2 UP, 3 LEFT)
    scores = [0, diag_score, up_score, left_score]
    best_score = max(scores)   # the 0 keeps scores from going negative, that's what makes it local
    # picks the first max on ties, so diag wins over gaps,
    # and END wins over diag when the diagonal score is exactly 0 (that's correct: the alignment stops there)
    direction  = np.argmax(scores)

    return best_score, direction


def traceback(seq1, seq2, traceback_matrix, maximum_position):
    """Rebuilds the local alignment by walking back from the best cell.

    Args:
        seq1 (str): sequence along the rows
        seq2 (str): sequence along the columns
        traceback_matrix (np.ndarray): direction codes from cal_score
        maximum_position (tuple): (row, col) of the highest score

    Returns:
        tuple: (aligned_seq1, aligned_seq2), same length, '-' for gaps
    """
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = maximum_position   # local alignment starts from the best scoring cell

    while traceback_matrix[i][j] != 0:   # 0 means END, the alignment starts here
        move = traceback_matrix[i][j]

        if move == 1:      # DIAG: both seqs give a letter
            aligned_seq1 += seq1[i-1]
            aligned_seq2 += seq2[j-1]
            i, j = i-1, j-1

        elif move == 2:    # UP: only the row moves, so seq2 gets the gap
            aligned_seq1 += seq1[i-1]
            aligned_seq2 += '-'
            i, j = i-1, j

        elif move == 3:    # LEFT: only the column moves, so seq1 gets the gap
            aligned_seq1 += '-'
            aligned_seq2 += seq2[j-1]
            i, j = i, j-1

    # we walked backwards from the end, so flip both strings
    return aligned_seq1[::-1], aligned_seq2[::-1]


def smith_waterman_align(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """Full week 5 Smith-Waterman: fill, find the best cell, traceback.

    Returns:
        tuple: (aligned_seq1, aligned_seq2, score_matrix)
    """
    # +1 for the padding row and column of zeros
    score_matrix     = np.zeros((len(seq1)+1, len(seq2)+1), dtype=int)
    traceback_matrix = np.zeros((len(seq1)+1, len(seq2)+1), dtype=int)

    max_score = 0
    max_pos = (0, 0)

    # start at 1 to skip the padding, each cell only needs neighbors we already filled
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score, direction = cal_score(score_matrix, seq1, seq2, i, j, match, mismatch, gap)
            score_matrix[i][j]     = score
            traceback_matrix[i][j] = direction

            # >= keeps the last cell that hits the top score if there's a tie
            if score >= max_score:
                max_score = score
                max_pos = (i, j)

    aligned_seq1, aligned_seq2 = traceback(seq1, seq2, traceback_matrix, max_pos)
    return aligned_seq1, aligned_seq2, score_matrix


if __name__ == "__main__":
    # test 1 from the week 5 slides (only runs when this file is run directly,
    # not when the notebook imports it)
    seq1 = 'TACTTAG'
    seq2 = 'CACATTAA'
    aligned_seq1, aligned_seq2, score_matrix = smith_waterman_align(seq1, seq2)
    print(aligned_seq1)
    print(aligned_seq2)
    print(score_matrix)
