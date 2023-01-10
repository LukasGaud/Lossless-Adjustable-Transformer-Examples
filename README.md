# Lossless-Adjustable-Transformer-Examples
## Installation
Copy the four .m files to the same location and simply run the scripts "thesisExampleQC.m" or "transformerModelling.m" in MATLAB. 

## thesisExampleQC.m
The script thesisExampleQC.m shows the calculations done in Example 5 of the thesis. 
Two functions fKron.m and fBlockDiag.m are needed to run the script (also included). They are summarised below:
  fKron.m - computes a transformation to convert the pencil sE-A into sT-K, such that both T and K are lower triangular
  fBlockDiag.m - takes the output from fKron.m and computes the Kronecker form of the regular matrix pencil sE-A
  
## transformerModelling.m
The script transformerModelling.m follows the methodology given in Section 3.5.2 and shows the steps when constructing a bilinear state-space representation for two RLC networks interconnected via lossless-adjustable transformer.

No additional functions are needed for this script.


