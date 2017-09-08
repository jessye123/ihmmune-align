package iHMMuneAlign;

/**
 * @author harris
 *
 * TODO
 * 
 */
/**
 * this class is used to calculate the mutability score of a penta nucleotide or
 * a trinucleotide. This class holds all the trinucleotide mutability scores
 */

public class MutabilityScoreTrinucleotide implements MutabilityScore {

	// Defines
	final char U = 'u';

	final char A = 'a';

	final char C = 'c';

	final char G = 'g';

	final char T = 't';

	// in order to look-up the nucleotide mutability scores quicker,
	// each nucleotide has a corresponding number which will be used
	// for direct index into an array (P1*16 + P2*4 + P3), where P1 is
	// the first nucleotide in the trinucleotide
	final int NUCL_A_NUMBER = 0;

	final int NUCL_C_NUMBER = 1;

	final int NUCL_G_NUMBER = 2;

	final int NUCL_T_NUMBER = 3;

	// array holding all the trinucleotide mutability score data
	// in order AAA,AAC,AAG...TTT
	/*
	 * final static double[] TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY =
	 * {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	 * 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	 * 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	 * 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	 * 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	 * 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0};
	 */
	// array holding all the adjusted trinucleotide mutability score data
	// in order AAA,AAC,AAG...TTT
	final static double[] TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY = { 0.65, 1.38,
			1.15, 2.01, 1.06, 1.14, 1.50, 1.29, 0.91, 2.12, 0.71, 1.32, 1.46,
			1.30, 1.18, 1.65, 1.17, 1.13, 1.51, 1.09, 1.08, 0.46, 0.75, 0.57,
			0.64, 0.56, 0.62, 1.25, 1.67, 0.60, 0.86, 0.64, 1.05, 0.49, 0.90,
			0.87, 1.73, 0.49, 0.30, 2.15, 0.71, 0.76, 0.71, 0.78, 3.06, 0.65,
			0.96, 1.29, 1.93, 1.76, 2.13, 1.61, 0.91, 0.69, 0.47, 0.61, 0.31,
			0.61, 0.58, 1.01, 1.32, 0.61, 0.28, 0.48 };

	// main method for testing only
	public static void main(String[] args) {
		MutabilityScoreTrinucleotide mutability_score_object = new MutabilityScoreTrinucleotide();

		System.out.println("TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY size = "
				+ TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY.length);

		double mutability_score = mutability_score_object
				.pentaNucleotideScore("cgtcc");

		System.out.println("mutability_score = " + mutability_score);
	}//--main()

	/**
	 * constructor
	 */
	public MutabilityScoreTrinucleotide() {

	}//--MutabilityScoreTrinucleotide

	/**
	 * calculates the mutability score of this pentanucleotide by taking the
	 * average mutability score of the three trinucleotides in this
	 * pentaNucleotide
	 * 
	 * returns the penta nucleotide mutatiblity score
	 */
	public double pentaNucleotideScore(String pentaNucleotide_any_case) {
		if (pentaNucleotide_any_case == null)
			throw new Error("penta nucleotide is null");

		// make the pentaNucleotide string lower case
		String pentaNucleotide = pentaNucleotide_any_case.toLowerCase();

		// pentanucleotide must have size 5
		if (pentaNucleotide.length() != 5)
			throw new Error(
					"pentaNucleotideScore(): pentaNucleotide is not of length 5");

		// total number of trinucleotides in one penta nucleotide
		final int TRINUCLEOTIDES_IN_PENTANUCLEOTIDE = 3;

		// chop this pentanucleotide into three trinucleotides,
		// find the nucleotide mutability score for each
		// and calculate the average (divide by three)

		String first_trinucleotide = pentaNucleotide.substring(0, 3);
		String middle_trinucleotide = pentaNucleotide.substring(1, 4);
		String last_trinucleotide = pentaNucleotide.substring(2, 5);

		// calculate mutability score for each trinucleotide
		double first_trinucleotide_score = triNucleotideScore(first_trinucleotide);

		double middle_trinucleotide_score = triNucleotideScore(middle_trinucleotide);

		double last_trinucleotide_score = triNucleotideScore(last_trinucleotide);

		// add up the trinucleotide scores
		double total_trinucleotide_score = first_trinucleotide_score
				+ middle_trinucleotide_score + last_trinucleotide_score;

		// get the average trinucleotide score
		double average_trinucleotide_score = total_trinucleotide_score
				/ TRINUCLEOTIDES_IN_PENTANUCLEOTIDE;

		// return the average trinucleotide score as a result
		return average_trinucleotide_score;

	}//--pentaNucleotideScore()

	/**
	 * finds the mutability score of tricucleotide "trinucleotide" the
	 * trinucleotide may have "U" symbols, which are unknown trinucleotides
	 * 
	 * returns mutability score of this trinucleotide
	 */
	private double triNucleotideScore(String trinucleotide) {
		// PPP System.out.println("triNucleotideScore()");

		// the resulting trinucleotide score
		double trinucleotide_score = Double.MIN_VALUE;

		// test the type of trinucleotide
		if (isNNN(trinucleotide)) {
			// PPP System.out.println("isNNN");

			// this is a normal trinucleotide, so we can lookup the score
			// directly
			trinucleotide_score = lookupTriNucleotideMutabilityScore(trinucleotide);
		} else if (isNNU(trinucleotide)) {
			// PPP System.out.println("isNNU");

			// the trinucleotide score is combination of scores
			trinucleotide_score = NNU_TriNucleotideScore(trinucleotide);
		} else if (isNUU(trinucleotide)) {
			// PPP System.out.println("isNUU");

			// the trinucleotide score is combination of scores
			trinucleotide_score = NUU_TriNucleotideScore(trinucleotide);
		} else if (isUNN(trinucleotide)) {
			// PPP System.out.println("isUNN");

			// the trinucleotide score is combination of scores
			trinucleotide_score = UNN_TriNucleotideScore(trinucleotide);
		} else if (isUUN(trinucleotide)) {
			// PPP System.out.println("isUUN");

			// the trinucleotide score is combination of scores
			trinucleotide_score = UUN_TriNucleotideScore(trinucleotide);
		}

		// return the resulting trinucleotide score
		return trinucleotide_score;

	}//--triNucleotideScore()

	/**
	 * looks-up the mutability score for trinucleotide "strict_trinucleotide",
	 * 
	 * Pre: "strict_trinucloeotide must contain a,c,g,t only, no "U" (unknown)
	 * nucleotides are allowed
	 * 
	 * returns mutability score of this trinucleotide
	 * 
	 * Erros if any of the nucleotides in "strict_trinucleotide" is not a,c,g or
	 * t
	 */
	private double lookupTriNucleotideMutabilityScore(
			String strict_trinucleotide) {
		// convert each trinucleotide in the "strict_trinucleotide" into
		// its corresponding number

		char charOne = strict_trinucleotide.charAt(0);
		int numberOne = getNucleotideNumber(charOne);

		char charTwo = strict_trinucleotide.charAt(1);
		int numberTwo = getNucleotideNumber(charTwo);

		char charThree = strict_trinucleotide.charAt(2);
		int numberThree = getNucleotideNumber(charThree);

		// now calculate the index into the TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY
		int trinucleotide_mutatbility_score_array_index = numberOne * 16
				+ numberTwo * 4 + numberThree;

		// lookup the mutability score using the calculated index
		double mutability_score = TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY[trinucleotide_mutatbility_score_array_index];

		// return the mutability score
		return mutability_score;

	}//--lookupTriNucleotideMutabilityScore()

	/**
	 * find the number corresponding to a nucleotide and return it
	 */
	private int getNucleotideNumber(char nucleotide) {
		switch (nucleotide) {
		case A:
			return NUCL_A_NUMBER;
		case C:
			return NUCL_C_NUMBER;
		case G:
			return NUCL_G_NUMBER;
		case T:
			return NUCL_T_NUMBER;
		default:
			throw new Error(
					"NucleotideNumber: nucleotide is not of type a,c,g ot t");
		}//--switch
	}//--getNucleotideNumber()

	/**
	 * find the nucleotide corresponding to a number and return it
	 */
	private char getNucleotideFromNumber(int nucleotide_number) {
		switch (nucleotide_number) {
		case NUCL_A_NUMBER:
			return A;
		case NUCL_C_NUMBER:
			return C;
		case NUCL_G_NUMBER:
			return G;
		case NUCL_T_NUMBER:
			return T;
		default:
			throw new Error(
					"getNucleotideFromNumber: nucleotide number is not in range");
		}//--switch
	}//--getNucleotideFromNumber()

	/* ********* Methods for testing trinucleotide type ********** */

	/**
	 * test if trinucleotide "trinucleotide" is of type NNN, which is a basic
	 * trinucleotide without any cut-offs
	 */
	private boolean isNNN(String trinucleotide) {
		// test each of the three characters in the trinucleotide separately

		char charOne = trinucleotide.charAt(0);
		if (!isN(charOne))
			return false;

		char charTwo = trinucleotide.charAt(1);
		if (!isN(charTwo))
			return false;

		char charThree = trinucleotide.charAt(2);
		if (!isN(charThree))
			return false;

		// all three nucleotides (characters) in trinucleotide had the correct
		// type, so return true
		return true;

	}//--isNNN()

	/**
	 * test if trinucleotide "trinucleotide" is of type NUU
	 */
	private boolean isNNU(String trinucleotide) {
		// test each of the three characters in the trinucleotide separately

		char charOne = trinucleotide.charAt(0);
		if (!isN(charOne))
			return false;

		char charTwo = trinucleotide.charAt(1);
		if (!isN(charTwo))
			return false;

		char charThree = trinucleotide.charAt(2);
		if (!isU(charThree))
			return false;

		// all three nucleotides (characters) in trinucleotide had the correct
		// type, so return true
		return true;
	}//--isNNU()

	/**
	 * test if trinucleotide "trinucleotide" is of type NUU
	 */
	private boolean isNUU(String trinucleotide) {
		// test each of the three characters in the trinucleotide separately

		char charOne = trinucleotide.charAt(0);
		if (!isN(charOne))
			return false;

		char charTwo = trinucleotide.charAt(1);
		if (!isU(charTwo))
			return false;

		char charThree = trinucleotide.charAt(2);
		if (!isU(charThree))
			return false;

		// all three nucleotides (characters) in trinucleotide had the correct
		// type, so return true
		return true;

	}//--isNUU()

	/**
	 * test if trinucleotide "trinucleotide" is of type NUU
	 */
	private boolean isUNN(String trinucleotide) {
		// test each of the three characters in the trinucleotide separately

		char charOne = trinucleotide.charAt(0);
		if (!isU(charOne))
			return false;

		char charTwo = trinucleotide.charAt(1);
		if (!isN(charTwo))
			return false;

		char charThree = trinucleotide.charAt(2);
		if (!isN(charThree))
			return false;

		// all three nucleotides (characters) in trinucleotide had the correct
		// type, so return true
		return true;

	}//--isUNN()

	/**
	 * test if trinucleotide "trinucleotide" is of type NUU
	 */
	private boolean isUUN(String trinucleotide) {
		// test each of the three characters in the trinucleotide separately

		char charOne = trinucleotide.charAt(0);
		if (!isU(charOne))
			return false;

		char charTwo = trinucleotide.charAt(1);
		if (!isU(charTwo))
			return false;

		char charThree = trinucleotide.charAt(2);
		if (!isN(charThree))
			return false;

		// all three nucleotides (characters) in trinucleotide had the correct
		// type, so return true
		return true;

	}//--isUUN()

	/* Helper Methods */

	/**
	 * test if nucleotide "nucleotide" is of type N (a,c,g or t)
	 * 
	 * returns true or false
	 */
	private boolean isN(char trinucleotide) {
		if (trinucleotide == A || trinucleotide == C || trinucleotide == G
				|| trinucleotide == T)
			return true;
		else
			return false;
	}//--isN

	/**
	 * test if nucleotide "nucleotide" is of type U ('u')
	 * 
	 * returns true or false
	 */
	private boolean isU(char trinucleotide) {
		if (trinucleotide == U)
			return true;
		else
			return false;
	}//--isU

	/* ********** End *********** */

	/**
	 * calculate the trinucleotide score of "NNU_trinucleotide" by taking the
	 * average of all combinations of trinucleotides possible to create, from
	 * this NNU type trinucleotide, by replacing the U symbol with a,c,g or t
	 * nucleotides
	 */
	private double NNU_TriNucleotideScore(String NNU_trinucleotide) {
		String temp_trinucleotide = ""; // current trinucleotide made by
										// replacement of U nucleotides
		char replacement_nucl; // the nucleotide replacement for the U symcol

		// sum of adding all the trinucleotide scores
		double total_trinucleotide_score = 0.0;

		// total number of possible trinucleotides that can be formed from
		// trinucleotide NNU
		final double TOTAL_TRINUCLEOTIDE_COMBINATIONS = 4;

		// make the 4 possible trinucleotide combinations by replacing the U at
		// the last position. combining/replacing from a..t
		for (int nucl_number = NUCL_A_NUMBER; nucl_number <= NUCL_T_NUMBER; nucl_number++) {
			// get the nucleotide corresponding to the nucleotide number
			replacement_nucl = getNucleotideFromNumber(nucl_number);

			// replace the U nucleotide with a,c,g,t nucleotide
			temp_trinucleotide = "" + NNU_trinucleotide.charAt(0)
					+ NNU_trinucleotide.charAt(1) + replacement_nucl;

			// PPP System.out.println("NNU_TriNucleotideScore produced: " +
			// temp_trinucleotide);

			// lookup the trinucleotide score for this NNN trinucleotide
			// combination
			double temp_trinucleotide_score = lookupTriNucleotideMutabilityScore(temp_trinucleotide);

			// add the trinucleotide combination score to the total score
			total_trinucleotide_score += temp_trinucleotide_score;

		}//--for(nucl_number)

		// find the average trinucleotide mutability score
		double average_trinucleotide_score = total_trinucleotide_score
				/ TOTAL_TRINUCLEOTIDE_COMBINATIONS;

		// return the average mutatbility score
		return average_trinucleotide_score;

	}//--NNU_TriNucleotideScore()

	/**
	 * calculate the trinucleotide score of "NUU_trinucleotide" by taking the
	 * average of all combinations of trinucleotides possible to create from
	 * this NNU type trinucleotide, by replacing the U symbol(s) with a,c,g or t
	 * nucleotides
	 */
	private double NUU_TriNucleotideScore(String NUU_trinucleotide) {
		// replace the middle nucleotide, so we get NNU, and then
		// call NNU for each of the four replacement combinations
		// and take the average of these four scores

		// number trinucleotide combinations which can be made from replacing
		// the
		// middle U in NUU (so NXU)
		final int TOTAL_NXU_TRINUCLEOTIDE_COMBINATIONS = 4;

		// the sum of all four NXU trinucleotide combination scores
		double total_NXU_trinucleotide_score = 0.0;

		// symbol X simply means the X will be replaced by a,c,g or t and thus
		// become N
		String temp_NXU_trinucleotide = ""; // current trinucleotide made by
											// replacement of middle U
											// nucleotides
		char middle_U_replacement_nucl; // the nucleotide replacement for the U
										// symcol in the mid position of the
										// trinucleotide

		// sum of adding all the trinucleotide scores
		double total_trinucleotide_score = 0.0;

		// make the 4 possible trinucleotide combinations by replacing the
		// middle U at
		// combining/replacing from a..t
		for (int nucl_number = NUCL_A_NUMBER; nucl_number <= NUCL_T_NUMBER; nucl_number++) {
			// get the nucleotide corresponding to the nucleotide number
			middle_U_replacement_nucl = getNucleotideFromNumber(nucl_number);

			// replace the middle U nucleotide with a,c,g,t nucleotide
			temp_NXU_trinucleotide = "" + NUU_trinucleotide.charAt(0)
					+ middle_U_replacement_nucl + NUU_trinucleotide.charAt(2);

			// PPP System.out.println("NUU_TriNucleotideScore produced: " +
			// temp_NXU_trinucleotide);

			// lookup the trinucleotide score for this NNU trinucleotide
			// combination
			double temp_NXU_trinucleotide_score = NNU_TriNucleotideScore(temp_NXU_trinucleotide);

			// add the trinucleotide NXU combination score to the total score
			total_NXU_trinucleotide_score += temp_NXU_trinucleotide_score;

		}//--for(nucl_number)

		// find the average trinucleotide mutability score
		double average_trinucleotide_score = total_NXU_trinucleotide_score
				/ TOTAL_NXU_TRINUCLEOTIDE_COMBINATIONS;

		// return the average mutatbility score
		return average_trinucleotide_score;

	}//--NUU_TriNucleotideScore()

	/**
	 * calculate the trinucleotide score of "NNU_trinucleotide" by taking the
	 * average of all combinations of trinucleotides possible to create from
	 * this NNU type trinucleotide, by replacing the U symbol with a,c,g or t
	 * nucleotides
	 */
	private double UNN_TriNucleotideScore(String UNN_trinucleotide) {
		String temp_trinucleotide = ""; // current trinucleotide made by
										// replacement of U nucleotides
		char replacement_nucl; // the nucleotide replacement for the U symcol

		// sum of adding all the trinucleotide scores
		double total_trinucleotide_score = 0.0;

		// total number of possible trinucleotides that can be formed from
		// trinucleotide NNU
		final double TOTAL_TRINUCLEOTIDE_COMBINATIONS = 4;

		// make the 4 possible trinucleotide combinations by replacing the U at
		// the first position. combining/replacing from a..t
		for (int nucl_number = NUCL_A_NUMBER; nucl_number <= NUCL_T_NUMBER; nucl_number++) {
			// get the nucleotide corresponding to the nucleotide number
			replacement_nucl = getNucleotideFromNumber(nucl_number);

			// replace the U nucleotide with a,c,g,t nucleotide
			temp_trinucleotide = "" + replacement_nucl
					+ UNN_trinucleotide.charAt(1) + UNN_trinucleotide.charAt(2);

			// PPP System.out.println("UNN_TriNucleotideScore produced: " +
			// temp_trinucleotide);

			// lookup the trinucleotide score for this NNN trinucleotide
			// combination
			double temp_trinucleotide_score = lookupTriNucleotideMutabilityScore(temp_trinucleotide);

			// add the trinucleotide combination score to the total score
			total_trinucleotide_score += temp_trinucleotide_score;

		}//--for(nucl_number)

		// find the average trinucleotide mutability score
		double average_trinucleotide_score = total_trinucleotide_score
				/ TOTAL_TRINUCLEOTIDE_COMBINATIONS;

		// return the average mutatbility score
		return average_trinucleotide_score;
	}//--UNN_TriNucleotideScore()

	/**
	 * calculate the trinucleotide score of "UUN_trinucleotide" by taking the
	 * average of all combinations of trinucleotides possible to create from
	 * this UUN type trinucleotide
	 */
	private double UUN_TriNucleotideScore(String UUN_trinucleotide) {
		// replace the middle nucleotide, so we get UNN, and then
		// call UNN for each of the four replacement combinations
		// and take the average of these four scores

		// number trinucleotide combinations which can be made from replacing
		// the
		// middle U in UUN (so UXN)
		final int TOTAL_UXN_TRINUCLEOTIDE_COMBINATIONS = 4;

		// the sum of all four UXN trinucleotide combination scores
		double total_UXN_trinucleotide_score = 0.0;

		// symbol X simply means the X will be replaced by a,c,g or t and thus
		// become N
		String temp_UXN_trinucleotide = ""; // current trinucleotide made by
											// replacement of middle U
											// nucleotides
		char middle_U_replacement_nucl; // the nucleotide replacement for the U
										// symcol in the mid position of the
										// trinucleotide

		// sum of adding all the trinucleotide scores
		double total_trinucleotide_score = 0.0;

		// make the 4 possible trinucleotide combinations by replacing the
		// middle U at
		// combining/replacing from a..t
		for (int nucl_number = NUCL_A_NUMBER; nucl_number <= NUCL_T_NUMBER; nucl_number++) {
			// get the nucleotide corresponding to the nucleotide number
			middle_U_replacement_nucl = getNucleotideFromNumber(nucl_number);

			// replace the middle U nucleotide with a,c,g,t nucleotide
			temp_UXN_trinucleotide = "" + UUN_trinucleotide.charAt(0)
					+ middle_U_replacement_nucl + UUN_trinucleotide.charAt(2);

			// PPP System.out.println("UUN_TriNucleotideScore produced: " +
			// temp_UXN_trinucleotide);

			// lookup the trinucleotide score for this UNN trinucleotide
			// combination
			double temp_UXN_trinucleotide_score = UNN_TriNucleotideScore(temp_UXN_trinucleotide);

			// add the trinucleotide UXN combination score to the total score
			total_UXN_trinucleotide_score += temp_UXN_trinucleotide_score;

		}//--for(nucl_number)

		// find the average trinucleotide mutability score
		double average_trinucleotide_score = total_UXN_trinucleotide_score
				/ TOTAL_UXN_TRINUCLEOTIDE_COMBINATIONS;

		// return the average mutatbility score
		return average_trinucleotide_score;
	}//--UUN_TriNucleotideScore()

}//--class MutabilityScoreTrinucleotide
