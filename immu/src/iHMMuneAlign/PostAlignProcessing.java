package iHMMuneAlign;
import java.util.Vector;
/**
 * rewrite by Jessica Ye
 * @author harris 
 */

/**
 * This class processes alignments after they are completed. The main purpose is
 * to filter out unsuccesfull alignments, and Mark them accordingly
 *  
 */
public class PostAlignProcessing {

	/*
	 * make sure a DGene alignment meets minimum requirements of a proper
	 * alignment If not, then replace all the DGene nucletodies with 'X', to
	 * show the user that this alignment has no proper DGene alignment
	 */
	public boolean acceptDGene(GeneInfo dGene, int dGeneAlignmentAcceptanceType) 
	{
		// compare the D-gene string with the aligned part of the sequence
		// string
		String geneString = dGene.aligned_gene_string;
		String alignmentString = dGene.aligned_ums_string;
		int alignmentLength = geneString.length();
		int mismatches = getMismatches(geneString, alignmentString);
		// largest number of consecutive matches in the DGene alignment
		int consecutiveMatches = getConsecutiveMatches(geneString, alignmentString);
		
		if(dGeneAlignmentAcceptanceType == GlobalDefines.EIGHT_MER_DGENE_ACCEPTANCE)
		{
			if(alignmentLength > 12)
				return analyzeTwelvePlusAlignment(geneString, alignmentString);
			else return hasEightMerDGeneAcceptance(alignmentLength, mismatches, consecutiveMatches);
		}
		else if(dGeneAlignmentAcceptanceType == GlobalDefines.FIVE_MER_CONSECUTIVE_DGENE_ACCEPTANCE)
		{
			return hasFiveMerDGeneAcceptance(alignmentLength, mismatches, consecutiveMatches);
		}
		else return false;
	}//acceptDGene()

	/**
	 * Test if the gene alignment fullfil the criteria for a specific alignment
	 * @param alignmentLength
	 * @param mismatches
	 * @param consecutiveMathces
	 * @return	true if accepted, othwerise false
	 */
	private boolean hasEightMerDGeneAcceptance(int alignmentLength, int mismatches, int consecutiveMathces)
	{
		if (alignmentLength < 8)
			return false;
		if (alignmentLength < 10) {
			// must have no mismatches
			if (mismatches <= 0)
				return true;
			else
				return false;
		}
		if (alignmentLength < 12) {
			// one mismatch is accepted
			if (mismatches <= 1)
				return true;
			else
				return false;
		}
		if (alignmentLength == 12) {
			// two mismatches is accepted
			if (mismatches <= 2)
				return true;
			else
				return false;
		}
		else return false;
	}//hasEightMerDGeneAcceptance()
	
	/**
	 * Test if the gene alignment fullfil the criteria for a specific alignment
	 * @param alignmentLength
	 * @param mismatches
	 * @param consecutiveMathces
	 * @return	true if accepted, otherwise false
	 */
	private boolean hasFiveMerDGeneAcceptance(int alignmentLength, int mismatches, int consecutiveMathces)
	{
		if (alignmentLength < 5)
			return false;
		else if(consecutiveMathces >= 5)
		{
			return true;
		}
		else return false;
	}//hasFiveMerDGeneAcceptance()
	
	/**
	 * For a gene string "geneString", aligned to "alignmentString", find number the largest
	 * number of consecutive matches
	 * @param geneString
	 * @param alignmentString
	 * @return		largest number of conescutive matches found
	 */
	private int getConsecutiveMatches(String geneString, String alignmentString)
	{
		char currGeneNucl;
		char currAlignedNucl;
		
		int consecutiveMatches = 0;
		int longestConsecutiveMatches = 0;
		
		// for every nucleotide in the gene string
		for(int i=0; i<geneString.length(); i++)
		{
			currGeneNucl = geneString.charAt(i);
			currAlignedNucl = alignmentString.charAt(i);
			// compare the two
			if(currGeneNucl == currAlignedNucl)
			{
				// match
				++consecutiveMatches;
			}
			else
			{
				// mismatch
				if(consecutiveMatches > longestConsecutiveMatches)
					longestConsecutiveMatches = consecutiveMatches;
				
				// reset the consecutive matches count
				consecutiveMatches = 0;
			}//else
		}//for
		
		// the string of consecutive matches may end in a match, so test for longest consecutive match
		if(consecutiveMatches > longestConsecutiveMatches)
			longestConsecutiveMatches = consecutiveMatches;
		
		return longestConsecutiveMatches;
	}//getConsecutiveMatches
	
	/*
	 * Look for the minimum requirements in a DGene EIGHT_MER sequence alignment
	 * where there are more than twelve nucloetides
	 */
	private boolean analyzeTwelvePlusAlignment(String geneString, String alignmentString)
	{
		int matches, mismatches;
		char geneNucl, alignmentNucl;

		// look for consecutive matches
		for (int i = 0; i < geneString.length(); i++) {
			matches = 0;
			mismatches = 0;
			for (int j = i; j < geneString.length(); j++) {
				geneNucl = geneString.charAt(i);
				alignmentNucl = alignmentString.charAt(i);

				// count mismatches and matches
				if (geneNucl != alignmentNucl) {
					++mismatches;
					if (mismatches > 2)
						break;
				}// nucl. mismatch
				else {
					++matches;
					if (matches >= 8) {
						if (mismatches == 0)
							return true;

						if (matches >= 10) {
							if (mismatches <= 1)
								return true;

							if (matches >= 12) {
								if (mismatches <= 2)
									return true;
							}// 12 or more matches
						}// 10 or more matches
					}// 8 or more matches

				}// nucleotide match

			}
		}

		return false; // min req. not found in this sequence
	}//analyzeTWELVE_PLUSSequence

	/**
	 * Get the total number of mismatches in a sequence
	 * @param geneString
	 * @param alignmentString
	 * @return
	 */
	private int getMismatches(String geneString, String alignmentString) {
		int mismatches = 0;
		char geneNucl, alignmentNucl;
		for (int i = 0; i < geneString.length(); i++) {
			geneNucl = geneString.charAt(i);
			alignmentNucl = alignmentString.charAt(i);

			if (geneNucl != alignmentNucl)
				mismatches++;
		}

		return mismatches;
	}//getMismatches
}