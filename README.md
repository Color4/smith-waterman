# Smith Waterman Algorithm

## Use


## Smith-Waterman implementation

### Picking best gap penalty and gap extension
After running through all the possible combinations of gap opening and gap extension penalties, the best false positive rate was 0.16, with a gap opening penalty of 4 and extension penalty of 4. I found this really interesting because I was expecting an opening penalty to be much larger than the extension penalty.

### Picking best scoring matrix 
At a true positive rate of 0.7, BLOSUM50 performed best. This makes sense to me as the gap penalties were optimized for BLOSUM50. I found that MATIO performed very poorly, which was expected. It was unexpected, however, that with these parameters, pam250 and blosum62 performed equally as poorly as MATIO. This suggests to me that perhaps the gap penalty of 4/4 was not a good indicator of alignment for all matrices.
