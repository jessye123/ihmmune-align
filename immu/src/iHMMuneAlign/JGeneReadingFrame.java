/*
 * Created on Aug 31, 2005
 *
 * TODO 
 * 
 */
package iHMMuneAlign;

/**
 * @author harris
 *
 * Methods to test if the JGene in a sequence is being read in the correct Frame
 */
public class JGeneReadingFrame {

	private boolean jGeneInFrame;
	private int WGXGmotifInAlignedJGenePos;

	// constructor
	public JGeneReadingFrame(String sequence, GeneInfo jGeneInfo, GeneInfo vGeneInfo) {
		
		jGeneInFrame = isJGeneInFrame(sequence, jGeneInfo, vGeneInfo);
		
	}//constructor
	
	/**
	 * 
	 * @param sequence
	 * @param jGeneInfo
	 * @param vGeneInfo
	 * @return
	 */
	private boolean isJGeneInFrame(String sequence, GeneInfo jGeneInfo, GeneInfo vGeneInfo)
	{
		// first, find the position in the aligned part of the JGene where
		// the WGXG pattern is located
		String alignedJGeneString = jGeneInfo.aligned_gene_string;
		int alignedJGeneWGXGmotifPos = findAlignedJGeneWGXGpos(alignedJGeneString); // counting from 0
		
		System.out.println("Aligned JGene motif position: " + alignedJGeneWGXGmotifPos);
		
		// now find the location of the aligned JGene in the sequence string,
		// and determine the exact postion of the WGXG pattern in the sequence string
		int jGeneStartPos = jGeneInfo.start_gene_pos; // counting from one
		
		// get the nucl. starting pos (in "sequence") of the part of the sequence which is aligned with the JGene
		int alignedJGeneSequencePos = sequence.indexOf(jGeneInfo.aligned_ums_string); 
		
		System.out.println("Aligned alignedJGeneSequencePos: " + alignedJGeneSequencePos);
		
		// add up the relative position of the JGene, the aligned JGene and WGXG motif position in the aligned JGene
		// to find the absolute position of the WGXG motif in the sequence
	//	int WGXGmotifIndexInSequence = ((jGeneStartPos - 1) + alignedJGeneSequencePos + alignedJGeneWGXGmotifPos);
		int WGXGmotifIndexInSequence = alignedJGeneSequencePos + alignedJGeneWGXGmotifPos;
		
		// now determine if the WGXG codon will be read in the correct frame based on the
		// reading frame of sequence's first codon (0,1 or 2)
		
		// determine the number of removed nucl. from the beginning of the VGene in the sequence
		int vGeneOffset = vGeneInfo.start_gene_pos - 1; // counting from 1
		// now calculate total number of nuceotides in preceding the WGXG motif from start of VGene
		int nuclPrecedingWGXGmotif = vGeneOffset + WGXGmotifIndexInSequence;
		
		System.out.println("Found the WGXG-motif at position: " + alignedJGeneWGXGmotifPos);
		this.WGXGmotifInAlignedJGenePos = alignedJGeneWGXGmotifPos;
		
		// if the number of nucleotides preceding the WGXG motif is can be divided by three,
		// (the size of a codon) then the JGene is read in-frame
		if(nuclPrecedingWGXGmotif % CodonHelperFunctions.CODON_LENGTH == 0)
		{
			// return the position of the WGXG motif in the aligned JGene part
			System.out.println("JGene in Frame");
			return true; 
		}
		// 
		System.out.println("JGene NOT in Frame");
		return false;
	}//f()
	
	/**
	 * Find the offset into the aligned part of the JGene where
	 * the WGXG pattern is located
	 * @param alignedJGeneString		nucleotide sequence string
	 * @return	the position of WGXG pattern in "alignedJGeneString"
	 * 			or -1 if not found
	 */
	private int findAlignedJGeneWGXGpos(String alignedJGeneString)
	{
		char[] alignedJGeneStringC = alignedJGeneString.toCharArray();
		char currNucl;
		char[] currCodon;
		// get/read every possible codon
		for(int i=0; i<alignedJGeneStringC.length - 2; i++)
		{
			currCodon = new char[]{alignedJGeneStringC[i], alignedJGeneStringC[i+1], alignedJGeneStringC[i+2]};
			if(CodonHelperFunctions.isAminoAcidW(currCodon))
			{
				System.out.println("found the W-codon in JGene at pos: " + i);
				// this codon is possibly the W codon in the WGXG motif, so look for the GXG codons
				if(hasGXGaminoAcids(i + 3, alignedJGeneStringC))
				{
					// the WGXG motif has been found in the aligned part of the JGene at position "i"
					return i;
				}
			}
		}
		return -1; // not found
	}//findAlignedJGeneWGXGpos
	
	/**
	 * Test if there if the amino acids XGX are locaed in the the nuclteotide sequence
	 * array "alignedJGeneStringC", at position "startPos"
	 * @param startPos
	 * @param alignedJGeneStringC
	 * @return	true if the GXG pattern is found, otherwise false
	 */
	private boolean hasGXGaminoAcids(int startPos, char[] alignedJGeneStringC)
	{
		// first test if we're at the end of the JGene String (meaning not enough
		// nucleotides availble for the GXG codon pattern
		int alignedJGeneLength = alignedJGeneStringC.length;
		int nucleotidesRequired = 3 * CodonHelperFunctions.CODON_LENGTH;
		int remainingNucleotidesInJGene = alignedJGeneLength - startPos;
		if(remainingNucleotidesInJGene < nucleotidesRequired)
			return false;
		
		// enough nucleotides for the GXG pattern, so test for it
		// get the three codons
		char[] codonOneG = new char[] {alignedJGeneStringC[startPos], alignedJGeneStringC[startPos + 1], alignedJGeneStringC[startPos + 2]};
		char[] codonTwoX = new char[] {alignedJGeneStringC[startPos + 3], alignedJGeneStringC[startPos + 4], alignedJGeneStringC[startPos + 5]};
		char[] codonThreeG = new char[] {alignedJGeneStringC[startPos + 6], alignedJGeneStringC[startPos + 7], alignedJGeneStringC[startPos + 8]};
		
		if(CodonHelperFunctions.isAminoAcidG(codonOneG)) // G
		{
			if(CodonHelperFunctions.isAminoAcidX(codonTwoX)) // X
			{
				if(CodonHelperFunctions.isAminoAcidG(codonThreeG)) // G
				{
					// all three codons match, so we have a successful match
					return true;
				}// if G
			}// if X
		}// if G
		
		return false; // no mathcing GXG pattern found
	}//hasGXGaminoAcids()

	public boolean isJGeneInFrame() {
		return jGeneInFrame;
	}
	public int getWGXGmotifInAlignedJGenePos() {
		return WGXGmotifInAlignedJGenePos;
	}
}//class JGeneReadingFrame
