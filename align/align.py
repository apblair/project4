# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub
    
    def _initialize_NeedlemanWunsch(self):
        """
        Initialize M, gap A, and gap B matrices.
        """
        self._align_matrix[0,0] = 0
        self._gapA_matrix[:,0] = [self.gap_open + n*self.gap_extend for n in range(len(self._seqA) + 1)] # column
        self._gapB_matrix[0,:] = [self.gap_open + m*self.gap_extend for m in range(len(self._seqB) + 1)] # row
        # self._initialize_NeedlemanWunsch_verbose()
    
    def _initialize_NeedlemanWunsch_verbose(self):
        """
        """
        print('### NW Setup ###', '\n')

        print('alignment matrix: ', '\n', self._align_matrix)
        print('gap a matrix: ', '\n', self._gapA_matrix)
        print('gap b matrix: ', '\n', self._gapB_matrix)
        print('\n')
        print('back trace: ', '\n', self._back)
        print('back a trace: ', '\n', self._back_A)
        print('back b trace', '\n', self._back_B)
        print('\n')

        print('seqA: ', self._seqA)
        print('seqB: ', self._seqB)
        print('\n')

        print('gap open: ', self.gap_open)
        print('gap extend: ',self.gap_extend)
        print('\n')

        print('### Initialization ###', '\n')
        print('alignment matrix: ', '\n', self._align_matrix)
        print('gap a matrix: ', '\n', self._gapA_matrix)
        print('gap b matrix: ', '\n', self._gapB_matrix)
        print('\n')

        print('### Compute inner score matrix ###', '\n')

    
    def _construct_M(self, m, n):
        """

        Parameters
        ----------
        m
        n
        """
        maxtrix_alignment_list = [self._align_matrix[m-1, n-1],
                                        self._gapA_matrix[m-1, n-1],
                                        self._gapB_matrix[m-1, n-1]]
        max_value, max_index = max(maxtrix_alignment_list), np.argmax(maxtrix_alignment_list)
        self._align_matrix[m,n] = self.sub_dict[(self._seqA[m-1], self._seqB[n-1])] + max_value
        self._back[m,n] = max_index
        # self._construct_M_verbose(m, n, maxtrix_alignment_list, max_value, max_index)
    
    def _construct_M_verbose(self, m, n, maxtrix_alignment_list, max_value, max_index):
        """

        Parameters
        ----------
        m
        n
        maxtrix_alignment_list
        max_value
        max_index
        """
        print('row: ', m)
        print('column: ', n, '\n')      

        print('## Alignment matrix ##')         
        print('matrix alignment list: ', maxtrix_alignment_list)
        print('max value: ', max_value)
        print('max index: ', max_index)
        print('amino acids: ', (self._seqA[m-1], self._seqB[n-1]))
        print('substitution dict value: ', self.sub_dict[(self._seqA[m-1], self._seqB[n-1])])
        print('alignment matrix: ', '\n', self._align_matrix)
        print('backtrace matrix: ', '\n', self._back, '\n')

    
    def _construct_gap_A(self, m, n):
        """

        Parameters
        ----------
        m
        n
        """
        a_matrix_list = [self.gap_open + self.gap_extend + self._align_matrix[m, n-1], 
                                self.gap_extend + self._gapA_matrix[m, n-1], 
                                self.gap_open + self.gap_extend + self._gapB_matrix[m, n-1]]
        a_max_value, a_max_index = max(a_matrix_list), np.argmax(a_matrix_list)
        self._gapA_matrix[m, n] = a_max_value
        self._back_A[m, n] = a_max_index
        # self._construct_gap_A_verbose(a_matrix_list, a_max_value, a_max_index)
        
    def _construct_gap_A_verbose(self, a_matrix_list, a_max_value, a_max_index):
        """

        Parameters
        ----------
        a_matrix_list
        a_max_value
        a_max_index
        """
        print('## A matrix ##')         
        print('A matrix list: ', a_matrix_list)
        print('max A value: ', a_max_value)
        print('max A index: ', a_max_index)
        print('gap A matrix: ', '\n',self._gapA_matrix)
        print('A backtrace matrix: ', '\n', self._back_A, '\n')

    def _construct_gap_B(self, m, n):
        """

        Parameters
        ----------
        m
        n
        """
        b_matrix_list = [self.gap_open + self.gap_extend + self._align_matrix[m-1,n],
                                 self.gap_open + self.gap_extend + self._gapA_matrix[m-1,n],
                                 self.gap_extend + self._gapB_matrix[m-1, n]] 
        b_max_value, b_max_index = max(b_matrix_list), np.argmax(b_matrix_list)
        self._gapB_matrix[m, n] = b_max_value
        self._back_B[m, n] = b_max_index
        # self._construct_gap_B_verbose(b_matrix_list, b_max_value, b_max_index)
    
    def _construct_gap_B_verbose(self, b_matrix_list, b_max_value, b_max_index):
        """

        Parameters
        ----------
        b_matrix_list
        b_max_value
        b_max_index
        """
        print('## B matrix ##')         
        print('B matrix list: ', b_matrix_list)
        print('max B value: ', b_max_value)
        print('max B index: ', b_max_index)
        print('gap B matrix: ', '\n', self._gapB_matrix)
        print('B backtrace matrix: ', '\n', self._back_B)
        print('\n','\n')
    
    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        Needleman-Wunsch algorithm to align protein or nucleotide sequences.
        
        Parameters
        ----------
        seqA : str
            A character sequence of amino acids or nucleotides.
        seqB : str
            A character sequence of amino acids or nucleotides.

        Returns
        -------
         tuple :
            a tuple of length 3 that contains the alignment score, alignment of sequence A, and alignment of sequence B.

        TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)
        Also, write up a docstring for this function using the
        _read_sub_matrix as an example.
        Don't forget to comment your code!
        """

        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        self._initialize_NeedlemanWunsch()

        # Calculate inner values in the score matrix
        for m in range(1, self._align_matrix.shape[0]): # iterate over sequence A (row)
            for n in range(1, self._align_matrix.shape[1]): # iterate over sequence B (column)
                self._construct_M(m,n)
                self._construct_gap_A(m,n)
                self._construct_gap_B(m, n)
        return self._backtrace()
    
    def _backtrace_verbose(self, m_row, n_column, max_matrix_score, max_matrix_index):
        """
        Parameters
        ----------
        m_row
        n_column
        max_matrix_score
        max_matrix_index)
        """
        print('### Starting backtrace ###')

        print('seq a: ', self._seqA)
        print('seq b: ',self._seqB)
        print('\n')

        print('Number of rows x columns: ', m_row, n_column)
        print('Alignment matrix shape: ', self._align_matrix.shape)
        print('\n')
        
        print('Max alignment score: ', max_matrix_score)
        print('Max alignment matrix index [align, gap a, gap b]: ', max_matrix_index)
        print('\n')

        print('Alignment seq A: ', self.seqA_align)
        print('Alignment seq B: ', self.seqB_align)
        print('Alignment score: ', self.alignment_score)

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        Traverse the three pointer matrices to construct the alignment of sequence A and B.

        Returns
        -------
        tuple :
            a tuple of length 3 that contains the alignment score, alignment of sequence A, and alignment of sequence B.
        
        References
        ----------
        1. https://wilkelab.org/classes/SDS348/2019_spring/labs/lab13-solution.html
        
        TODO: Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        # Implement this method based upon the heuristic chosen in the align method above.

        m_row,n_column = len(self._seqA),len(self._seqB)

        matrices_list = [self._align_matrix[-1,-1], self._gapA_matrix[-1,-1], self._gapB_matrix[-1,-1]]
        max_matrix_score, max_matrix_index  = max(matrices_list), np.argmax(matrices_list)
        self.alignment_score = max_matrix_score

        while m_row>0 and n_column>0:
            
            if max_matrix_index == 0:
                self.seqA_align = self._seqA[m_row-1] + self.seqA_align 
                self.seqB_align = self._seqB[n_column-1] + self.seqB_align 
                max_matrix_index = self._back[m_row,n_column]
                m_row-=1
                n_column-=1

            elif max_matrix_index == 1:
                self.seqA_align = '-' + self.seqA_align 
                self.seqB_align = self._seqB[n_column-1] + self.seqB_align 
                max_matrix_index = self._back_A[m_row, n_column]
                n_column-=1

            elif max_matrix_index == 2: 
                self.seqA_align = self._seqA[m_row-1]+ self.seqA_align
                self.seqB_align = '-' + self.seqB_align
                max_matrix_index = self._back_B[m_row, n_column]
                m_row-=1
        
        # self._backtrace_verbose(m_row, n_column, max_matrix_score, max_matrix_index)

        return (self.alignment_score, self.seqA_align, self.seqB_align)

def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header