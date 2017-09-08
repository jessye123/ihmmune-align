package iHMMuneAlign;

/**
 * @author harris
 *
 * TODO
 * 
 */

/**
 * this class holds methods for calculating the mutability score of a
 * nucleotide, based on on whether it is part of a mutational hot spot or not.
 */
public class MutabilityScoreHotspot implements MutabilityScore {
	// /////////////////// Defines

	// //////// Nucleotides A,C,G,T
	final char A_NUCLEOTIDE = 'a';

	final char C_NUCLEOTIDE = 'c';

	final char G_NUCLEOTIDE = 'g';

	final char T_NUCLEOTIDE = 't';

	final char N_NUCLEOTIDE = 'n'; // in this context, a n nucleotide we are
								   // unable to give

	// an a,c,g,t nucleotide, becuase we don't know (gene end)

	// the probability of a mutation at a nt. at pos "p" is already calculated,
	// we only find the mutability score for the hot spots (observed mutations /
	// expected mutations)

	// Pre calculated for simplicity
	final double ONE_THIRD = ((double) 1 / (double) 3);

	final double ONE_SIXTH = ((double) 1 / (double) 6);

	final double ONE_EIGTH = ((double) 1 / (double) 8);

	final double ONE_QUARTER = ((double) 1 / (double) 4);

	final double THREE_QUARTER = ((double) 3 / (double) 4);

	final double ONE_OVER_32 = ONE_EIGTH * ONE_QUARTER;

	// Hotspot nucleotide coverage
	// - The percentage of nucleotides belonging to a specfic hotspot
	final double RGYW_Nucleotide_coverage = ONE_OVER_32;

	final double WRCY_Nucleotide_coverage = ONE_OVER_32;

	final double WAN_Nucleotide_coverage = ONE_EIGTH;

	// the percentage of nucleotide NOT belonging to any hotspots
	final double NO_HOTSPOT_Nucleotide_coverage = 1.00 - (RGYW_Nucleotide_coverage
			+ WRCY_Nucleotide_coverage + WAN_Nucleotide_coverage);

	// Mutation probabilities scores (percentage wise) for GYW, WRCY and WAN
	final double RGYW_mutation_prob = 0.138;

	final double WRCY_mutation_prob = 0.096;

	final double WAN_mutation_prob = 0.161;

	// the mutation score for nucleotides not part of any hotspots
	final double NO_HOTSPOT_mutation_prob = (1.00 - (RGYW_mutation_prob
			+ WRCY_mutation_prob + WAN_mutation_prob));

	// original mutation probability scores
	//final double RGYW_mutation_prob = ONE_THIRD;
	//final double WRCY_mutation_prob = 0.096;
	//final double WAN_mutation_prob = 0.161;

	// The mutability score for the Non Hot Spots is calculated just like a hot
	// spot
	final double NO_HOTSPOT_MUTABILITY_SCORE = NO_HOTSPOT_mutation_prob
			/ NO_HOTSPOT_Nucleotide_coverage;

	// the actual mutability score for a hot spot is calculated by taking the
	// nucleotide mutation score
	// and divide it by the nucleotide coverage score.
	final double WAN_MUTABILITY_SCORE = WAN_mutation_prob
			/ WAN_Nucleotide_coverage;

	// there are special cases when the quad nucleotide can not be determined
	// because
	// it contains N-nucleotide (it's at the start or end of a gene), for these
	// cases we pre calculate the mutability scores
	// 50% chance the NAN is a WAN
	final double NAN_MUTABILITY_SCORE = ((WAN_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE) / 2);

	final double RGYW_MUTABILITY_SCORE = RGYW_mutation_prob
			/ RGYW_Nucleotide_coverage;

	final double NGYW_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE) / 2);

	final double RGYN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE) / 2);

	final double RGNN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE * ONE_QUARTER) + (NO_HOTSPOT_MUTABILITY_SCORE * THREE_QUARTER));

	// WRCY covers (1/2 * 1/2 * 1/4 * 1/2) 1/32 of all quadnucleotides, but
	// covers 1/6 of all mutations, therefore, the mutability score for a "c"
	// nt.
	// which is part of a RGYW quadnucleotide is:
	//	final double ONE_OVER_NINE = ((double)1/(double)9);

	final double WRCY_MUTABILITY_SCORE = WRCY_mutation_prob
			/ WRCY_Nucleotide_coverage;

	final double NRCY_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE) / 2);

	final double NNCY_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE * ONE_QUARTER) + (NO_HOTSPOT_MUTABILITY_SCORE * THREE_QUARTER));

	final double WRCN_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE) / 2);

	// ////////////// Globals
	boolean DEBUGGING = false; // are we debugging or not

	/**
	 * constructor
	 */
	public MutabilityScoreHotspot() {

	}//--MutabilityScoreHotspot()

	/**
	 * get the mutability score for a nucleotide in the center of penta
	 * nucleotide "pentaNucleotide". "pentanNucleotide" can contain the
	 * NUCLEOTIDE 'n', indicating we are at the end of a gene (f.ex: nnagt,
	 * means that we are at the first position in a gene)
	 */
	public double pentaNucleotideScore(String pentaNucleotideUnknownCase) {
		// make the pentaNucleotide string into lower case, as all nts.
		// definitions in this
		// class is in lower case
		String pentaNucleotide = pentaNucleotideUnknownCase.toLowerCase();

		// get the center nucleotide from the penta nucleotide
		char centerNucleotide = pentaNucleotide.charAt(2); // 0,1,2,3,4 (2 is he
														   // middle of the
														   // string)

		// test if the nucleotide is part of any hotspots
		switch (centerNucleotide) {
		case T_NUCLEOTIDE: // no hotspots are centered around a T nucleotide
			if (DEBUGGING)
				System.out.println("Main: No Hot");
			return NO_HOTSPOT_MUTABILITY_SCORE; // so return the no hot spot
												// mutability score
		case A_NUCLEOTIDE:
			return testWAN(pentaNucleotide);
		case G_NUCLEOTIDE:
			return testRGYW(pentaNucleotide);
		case C_NUCLEOTIDE:
			return testWRCY(pentaNucleotide);
		}//--switch

		throw new Error(
				"pentaNucleotideScore(): center nucleotide not a,c,g or t, but : "
						+ centerNucleotide);

	}//--pentaNucleotideScore()

	/**
	 * test if this "a" nucleotide centered pentaNucleotide is actually part of
	 * the WAN trinucleotide returns WAN_MUTATBILITY_SCORE if above test is true
	 * else return NO_HOTSPOT_MUTABILITY_SCORE or if pentaNucleotide is part of
	 * a gene end, return approximation
	 */
	private double testWAN(String pentaNucleotide) {
		// get the second nucleotide (W) from the pentaNucleotide
		char secondNucleotide = pentaNucleotide.charAt(1);

		// if the second nucleotide is W, we have a WAN trinucleotide
		if (isW(secondNucleotide)) {
			if (DEBUGGING)
				System.out.println("WAN");
			return WAN_MUTABILITY_SCORE; // return the WAN mutability score as
										 // this is a WAN trinucleotide
		} else {
			// test if the second nucleotide is a N
			if (isN(secondNucleotide)) {
				// since there is a 50/50 probability that the N (second)
				// nucleotide would be
				// a W, return the combined mutability score of no-hotspot and
				// WAN
				if (DEBUGGING)
					System.out.println("NAN");
				return NAN_MUTABILITY_SCORE;
			}//--if
			else // there is not hotspot, and we are not at an end of a gene
			{
				// return the no-hotspot mutability score
				if (DEBUGGING)
					System.out.println("WAN: No Hot");
				return NO_HOTSPOT_MUTABILITY_SCORE;

			}//--e;se
		}//--else

	}//--testnWANn()

	/**
	 * test if this "g" nucleotide centered pentaNucleotide is actually part of
	 * the RGYW quadnucleotide returns RGYW_MUTATBILITY_SCORE if above test is
	 * true else return NO_HOTSPOT_MUTABILITY_SCORE or if pentaNucleotide is
	 * part of a gene end, return approximation
	 */
	private double testRGYW(String pentaNucleotide) {
		// get the second (R), fourth (Y) and fifth (W) nucleotide from the
		// pentaNucleotide
		char secondNucleotide = pentaNucleotide.charAt(1);
		char fourthNucleotide = pentaNucleotide.charAt(3);
		char fifthNucleotide = pentaNucleotide.charAt(4);

		// in order for this to be a RGYW quad nucleotide, second, fourth and
		// fifth nt. need
		// to be correct
		if (isR(secondNucleotide) && isY(fourthNucleotide)
				&& isW(fifthNucleotide)) {
			return RGYW_MUTABILITY_SCORE; // is RGYW quad nucleotide, so return
										  // the RGYW mutability score
		}//--if
		else //  this is not a RGYW quad nucleotide, so test if this is due to
			 // any n nucleotides
		{
			if (isN(secondNucleotide)) // is this the start end of a gene (with
									   // first nucleotide equal N)
			{
				if (isY(fourthNucleotide) && isW(fifthNucleotide)) {
					// we have a quad nucleotide of type NGYW
					// so return the pre-calculated NGYW mutability score
					if (DEBUGGING)
						System.out.println("NGYW");
					return NGYW_MUTABILITY_SCORE;

				}//--if

			}//--if
			else if (isN(fourthNucleotide)) // is this the end of a gene (with
											// last two nts. equal N)
			{
				if (isR(secondNucleotide)) {
					// we have a quad nucleotide of type RGNN
					// so return the pre-calculated RGNN mutability score
					if (DEBUGGING)
						System.out.println("RGNN");
					return RGNN_MUTABILITY_SCORE;

				}//--if

			}//--else if
			else if (isN(fifthNucleotide)) // is this the end of a gene
			{
				if (isR(secondNucleotide) && isY(fourthNucleotide)) {
					// we have a quad nucleotide of type RGYN
					// so return the pre-calculated RGYN mutability score
					if (DEBUGGING)
						System.out.println("RGYN");
					return RGYN_MUTABILITY_SCORE;

				}//--if

			}//--else if
		}//else

		// there is not hotspot, and we are not at an end of a gene
		// return the no-hotspot mutability score
		if (DEBUGGING)
			System.out.println("RGYW: No Hot");

		return NO_HOTSPOT_MUTABILITY_SCORE;

	}//--testRGYW()

	/**
	 * test if this "c" nucleotide centered pentaNucleotide is actually part of
	 * the WRCY quadnucleotide returns WRCY_MUTATBILITY_SCORE if above test is
	 * true else return NO_HOTSPOT_MUTABILITY_SCORE or if pentaNucleotide is
	 * part of a gene end, return approximation
	 */
	private double testWRCY(String pentaNucleotide) {
		// get the first (W), second (R) and fourth (Y) nucleotide from the
		// pentaNucleotide
		char firstNucleotide = pentaNucleotide.charAt(0);
		char secondNucleotide = pentaNucleotide.charAt(1);
		char fourthNucleotide = pentaNucleotide.charAt(3);

		// in order for this to be a RGYW quad nucleotide, first, second and
		// fourth nt. need
		// to be correct
		if (isW(firstNucleotide) && isR(secondNucleotide)
				&& isY(fourthNucleotide)) {
			if (DEBUGGING)
				System.out.println("WRCY");
			return WRCY_MUTABILITY_SCORE; // is WRCY quad nucleotide, so return
										  // the WRCY mutability score
		}//--if
		else //  this is not a WRCY quad nucleotide, so test if this is due to
			 // any n nucleotides
		{
			if (isN(fourthNucleotide)) // is this the end of a gene (last two
									   // are N)
			{
				if (isW(firstNucleotide) && isR(secondNucleotide)) {
					// we have a quad nucleotide of type WRCN
					// so return the pre-calculated WRCN mutability score
					if (DEBUGGING)
						System.out.println("WRCN");
					return WRCN_MUTABILITY_SCORE;

				}//--if

			}//--if
			else if (isN(secondNucleotide)) // is this the start end of a gene
											// (first two are N)
			{
				if (isY(fourthNucleotide)) {
					// we have a quad nucleotide of type NNCY
					// so return the pre-calculated NNCY mutability score
					if (DEBUGGING)
						System.out.println("NNCY");
					return NNCY_MUTABILITY_SCORE;

				}//--if

			}//--else if
			else if (isN(firstNucleotide)) // is this the start end of a gene
										   // (only first nucleotide is N)
			{
				if (isR(secondNucleotide) && isY(fourthNucleotide)) {
					// we have a quad nucleotide of type NRCY
					// so return the pre-calculated NRCY mutability score
					if (DEBUGGING)
						System.out.println("NRCY");
					return NRCY_MUTABILITY_SCORE;

				}//--if

			}//--else if
		}//else

		// there is not hotspot, and we are not at an end of a gene
		// return the no-hotspot mutability score
		if (DEBUGGING)
			System.out.println("WRCY: No Hot");
		return NO_HOTSPOT_MUTABILITY_SCORE;

	}//--testWRCY()

	/**
	 * test if nucleotide "nucleotide" is of type N, N_NUCLEOTIDE is only used
	 * at the end of a gene when a whole penta nucleotide is not available
	 * returns true or false
	 */
	private boolean isN(char nucleotide) {
		if (nucleotide == N_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isA

	/**
	 * test if nucleotide "nucleotide" is of type A returns true or false
	 */
	private boolean isA(char nucleotide) {
		if (nucleotide == A_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isA

	/**
	 * test if nucleotide "nucleotide" is of type C returns true or false
	 */
	private boolean isC(char nucleotide) {
		if (nucleotide == C_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isC

	/**
	 * test if nucleotide "nucleotide" is of type G returns true or false
	 */
	private boolean isG(char nucleotide) {
		if (nucleotide == G_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isG

	/**
	 * test if nucleotide "nucleotide" is of type T returns true or false
	 */
	private boolean isT(char nucleotide) {
		if (nucleotide == T_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isT

	// //////////////////////////////////////////////////////

	/**
	 * test if nucleotide "nucleotide" is of type R (A/G) returns true or false
	 */
	private boolean isR(char nucleotide) {
		if (isA(nucleotide) || isG(nucleotide))
			return true;
		else
			return false;
	}//--isR

	/**
	 * test if nucleotide "nucleotide" is of type Y (C/T) returns true or false
	 */
	private boolean isY(char nucleotide) {
		if (isC(nucleotide) || isT(nucleotide))
			return true;
		else
			return false;
	}//--isY

	/**
	 * test if nucleotide "nucleotide" is of type W (A/T) returns true or false
	 */
	private boolean isW(char nucleotide) {
		if (isA(nucleotide) || isT(nucleotide))
			return true;
		else
			return false;
	}//--isW

	/**
	 * for testing only
	 */
	public static void main(String[] args) {
		MutabilityScoreHotspot object = new MutabilityScoreHotspot();

		System.out.println("dividing two integers: 5 / 2 = " + (double) 5 / 2);

		double mutabilityScore = object.pentaNucleotideScore("naacg");
		  	mutabilityScore = object.pentaNucleotideScore("agccn");
			mutabilityScore = object.pentaNucleotideScore("ttctn");
		// 	mutabilityScore = object.pentaNucleotideScore("cagtt");
		// 	mutabilityScore = object.pentaNucleotideScore("nactn");
		// 	mutabilityScore = object.pentaNucleotideScore("acgtn");

		System.out.println("ONE_EIGTH mutability score = " + object.ONE_EIGTH);
		System.out.println("WAN mutability score = "
				+ object.WAN_MUTABILITY_SCORE);
		System.out.println("mutability score = " + mutabilityScore);
	}//--main()

}//--MutabilityScoreHotspot
