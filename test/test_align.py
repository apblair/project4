# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    alignment_tuple = nw.align(seq1, seq2)

    assert tuple([len(seq1)+1, len(seq2)+1]) == nw._align_matrix.shape, "Error: matrix dimensions are wrong. The shape should be rows of length of seq1+1 and columns of length of seq2+1."
    assert nw._align_matrix[0, 0] == 0, "Error: Initialization was not set correctly!"
    assert alignment_tuple[0] == 4, "Error: Alignment score for test case is wrong!"
    assert np.allclose(nw._align_matrix[-1,:],[-np.inf, -14, -6, 4])

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    pass


test_nw_alignment()

