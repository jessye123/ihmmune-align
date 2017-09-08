package iHMMuneAlign;


import java.io.PrintWriter;

/**
 * this class is a container for one or more static methods used to calculate
 * the A value for a given VGene The A-value for an alignment is the probability
 * of an alignment in the VGene at pos 0, this probability decays exponentially
 * as we move away from the start of the VGene
 */
public class A_Score {
	/**
	 * the A value used for the CDR prediction is found by looking at the number
	 * of mutations in the common region of a sequence over the number of
	 * nucleotides in this common region multiplied by the A_slope for this
	 * common region.
	 * 
	 * The A_slope is previously determined by using linear regression on the A
	 * values (Y's) over the mutation ratio for the common region (X's) and thus
	 * estimating the relationship between A value and number of of mutations
	 * for common region
	 */

	// the n/x nts which may be found in some sequences
	final static char N_NUCLEOTIDE = 'n';
	final static char X_NUCLEOTIDE = 'x';

	// the start and end position of the common region (starting at 1)
	final static int COMMON_AREA_START_POS = 129;
	final static int COMMON_AREA_END_POS = 243;
	final static int COMMON_AREA_NUCL_LENGTH = COMMON_AREA_END_POS
			- COMMON_AREA_START_POS + 1;

	// the A-slope value previously determined mathcing the
	// common region specified above
	final static double A_SLOPE = 0.0065;
	final static double B_ADDITION = 0.0015; // this is the b found from (ax +b)
											 

	/**
	 * calculate the A probability using the above rules and the number of
	 * mutations found in the V-Gene "VGene" when compared to it's mathcing
	 * sequence "matching sequence", the "VGene_start_offset" must be taken into
	 * account when placing the common region
	 * 
	 * Pre: "VGene" and "matching_sequence" are aligned
	 * 
	 * returns A probability based on number of mutations in VGene's common area
	 * plus above rules
	 * 
	 * Errors if VGene does not contain the entire common region, and if
	 * matching sequence is not at least the length of the VGene
	 */
	public static String saySomething (){
		return "Hello World";
	}
	public static double A_probability(String VGeneUnknownCase,
			String matching_sequenceUnknownCase, int VGene_start_offset,
			PrintWriter pw) {
		// get the lower case string version of the method variables
		String VGene = VGeneUnknownCase.toLowerCase();
		String matching_sequence = matching_sequenceUnknownCase.toLowerCase();

		// ////////// ERROR TESTING (ASSERTIONS)

		System.out.println("VGene length = " + VGene.length());
		System.out.println("Matching seq. length = "
				+ matching_sequenceUnknownCase.length());

		// does V-Gene contain common region (currently is must, but this can
		// easily be changed)

		if ((VGene.length() + VGene_start_offset) < COMMON_AREA_END_POS)
			throw new Error(
					"A_Score: A_probability(): VGene does not contain common region (VGene too short)");
		if (VGene_start_offset > COMMON_AREA_START_POS)
			throw new Error(
					"A_Score: A_probability(): VGene does not contain common region (VGene offset to great)");

		// VGene Length must be shorter than or equal to it's matching sequence
		// length
		if (VGene.length() > matching_sequence.length())
			throw new Error(
					"A_Score: A_probability(): VGene length is greater than its matching sequence length");

		/////////////////// calculation of probability of mutation using the
		// number of mutations in the VGene

		// find number of mutatons in VGene
		int mutationsInVGene = mutationsInVGene(VGeneUnknownCase,
				matching_sequenceUnknownCase, VGene_start_offset, pw);

		// find the length of the VGene (the start of the VGene is often
		// missing)
		int vGeneLength = VGeneUnknownCase.length();

		// get the ratio between of CR length over VGene length
		double CRtoVGeneLengthRatio = (double) COMMON_AREA_NUCL_LENGTH
				/ (double) vGeneLength;

		// multiply the number of mutations in the VGene with the length ratio,
		// and we have "averaged" number of mutations in the CR
		double averageMutationsInCR = (double) mutationsInVGene
				* (double) CRtoVGeneLengthRatio;

		// ///////////////////////// find number of mutations in VGene (by
		// looking at the matching sequence)

		int number_of_mutations_in_CR = 0; // CR = common region
		int number_of_NandX_nts_in_CR = 0; // the number of N/X nts in the
										   // common region, these are not
										   // actually mutations, but could be

		char V_nucl, MS_nucl; // current VGene and mathcing sequence Nucleotide

		// for every nucleotide in the Common Region of the VGene
		// "VGene_start_offset" must be taken into account
		// for common area start and end position
		for (int i = (COMMON_AREA_START_POS - 1 - VGene_start_offset); i <= (COMMON_AREA_END_POS - 1 - VGene_start_offset); i++) {
			V_nucl = VGene.charAt(i);
			MS_nucl = matching_sequence.charAt(i);

			// test if the nucleotide's are not the same
			if (V_nucl != MS_nucl) {
				// test if the difference is caused by a n/x nucleotide,
				// in which case it does not count as a full mutation
				if (MS_nucl == N_NUCLEOTIDE || MS_nucl == X_NUCLEOTIDE) {
					++number_of_NandX_nts_in_CR; // found a n/x nucleotide, so
												 // increment count
				} else {
					// a regular mutation occured
					++number_of_mutations_in_CR; // found a mutation, so
												 // increase count
				}//--else
			}//--if
		}//--for(i)

		System.out.println("number of mutations found in common area = "
				+ number_of_mutations_in_CR);
		System.out.println("number of n/x nts. found in common area = "
				+ number_of_NandX_nts_in_CR);

		// if the Print Writer is not null, write number of mutations found inCR
		// to file
		if (pw != null) {
			pw.println(number_of_mutations_in_CR + "\t"
					+ COMMON_AREA_NUCL_LENGTH);
		}//--if

		////////////////// The calcualation of the Probability of mutation at
		// any nucleotide
		////////////////// given the number of mutations in the CR

		// when calculating the A-probability, I wanna take into account the
		// number of n/x nts found as well,
		// so based on how many real mutations I found, the likelihood of a n/x
		// being a mutation increases
		// based on the formula: ((M * A_SLOPE) + B_ADDITION) + ((M / L) * U *
		// A_SLOPE)
		// where
		// M = "real mutations found"
		// L = "length of common area"
		// U = "number of N/X nts found"
		// B = "the b in AX + B"
		// we get the complete A_Probability

		// calculate A based on "mutations_in_common_region" and A_slope
		// B + (M * A_Slope)
		double probabilityMutationPlusB = ((number_of_mutations_in_CR * A_SLOPE) + B_ADDITION);

		// either the averaged (estimation) of number of mutations in the CR, or
		// the actual number of mutations in the CR
		// can be used to get the probability of a mutation in a nucleotide
		//	    double NandXaddition = (( (double)number_of_mutations_in_CR /
		// (double)COMMON_AREA_NUCL_LENGTH)* A_SLOPE *
		// number_of_NandX_nts_in_CR);
		double NandXaddition = (((double) averageMutationsInCR / (double) COMMON_AREA_NUCL_LENGTH)
				* A_SLOPE * number_of_NandX_nts_in_CR);

		double aProbability = probabilityMutationPlusB + NandXaddition;

		// return result
		return aProbability;

	}//--A_probability()

	/**
	 * get the number of mutations in the aligned part of the VGene (complete
	 * VGene minus start offset) Start offset has already been removed by the
	 * alignment process
	 */
	private static int mutationsInVGene(String VGeneUnknownCase,
			String matching_sequenceUnknownCase, int VGene_start_offset,
			PrintWriter pw) {
		// get the lower case string version of the method variables
		String VGene = VGeneUnknownCase.toLowerCase();
		String matching_sequence = matching_sequenceUnknownCase.toLowerCase();

		// find number of mutations in VGene (by looking at the matching
		// sequence)

		int number_of_mutations_in_VGene = 0; // CR = common region
		int number_of_NandX_nts_in_VGene = 0; // the number of N/X nts in the
											  // common region, these are not
											  // actually mutations, but could
											  // be

		char V_nucl, MS_nucl; // current VGene and mathcing sequence Nucleotide

		// for every nucleotide in the the VGene
		for (int i = 0; i < VGene.length(); i++) {
			V_nucl = VGene.charAt(i);
			MS_nucl = matching_sequence.charAt(i);

			// test if the nucleotide's are not the same
			if (V_nucl != MS_nucl) {
				// test if the difference is caused by a n/x nucleotide,
				// in which case it does not count as a full mutation
				if (MS_nucl == N_NUCLEOTIDE || MS_nucl == X_NUCLEOTIDE) {
					++number_of_NandX_nts_in_VGene; // found a n/x nucleotide,
													// so increment count
				} else {
					// a regular mutation occured
					++number_of_mutations_in_VGene; // found a mutation, so
													// increase count
				}//--else
			}//--if
		}//--for(i)

		int vGeneLength = VGene.length();

		System.out.println("number of mutations found in VGene area = "
				+ number_of_mutations_in_VGene);
		System.out.println("number of n/x nts. found in VGene  area = "
				+ number_of_NandX_nts_in_VGene);
		System.out.println("VGene area length = " + vGeneLength);

		// if the Print Writer is not null, write number of mutations found inCR
		// to file
		if (pw != null) {
			pw.print(number_of_mutations_in_VGene + "\t" + VGene_start_offset
					+ "\t" + vGeneLength + "\t");
		}//--if

		// return number of mutations found in the VGene
		return number_of_mutations_in_VGene;

	}//--mutationsInVGene()

	/**
	 * main method for testing
	 */
	public static void main(String[] args) {
		System.out.println("test begin");
		double A_prob = A_probability(
				"aaacccgggtttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaagaatacaa",
				"aaacccaaatttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaggggggggggggggggggggggggggggggg",
				128, null);

		System.out.println("A_prob = " + A_prob);
	}//--main

}//--A_Score
