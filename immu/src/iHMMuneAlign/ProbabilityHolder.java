package iHMMuneAlign;

/**
 * @author harris
 * 
 * TODO
 *  
 */
import java.util.ArrayList;
import java.util.StringTokenizer;

public class ProbabilityHolder {
	private static final String VD_N = "VD_N";

	private static final String DJ_N = "DJ_N";

	private static final String V_END_EXO_MEAN = "V_end_exo_mean";

	private static final String V_END_EXO_STDDEV = "V_end_exo_stdDev";

	private static final String D_START_EXO_MEAN = "D_start_exo_mean";

	private static final String D_START_EXO_STDDEV = "D_start_exo_stdDev";

	private static final String D_END_EXO_MEAN = "D_end_exo_mean";

	private static final String D_END_EXO_STDDEV = "D_end_exo_stdDev";

	private static final String J_START_EXO_MEAN = "J_start_exo_mean";

	private static final String J_START_EXO_STDDEV = "J_start_exo_stdDev";

	private static final String V_END_P = "V_end_P";

	private static final String D_START_P = "D_start_P";

	private static final String D_END_P = "D_end_P";

	private static final String J_START_P = "J_start_P";

	private static final String GENE_MUTATION = "Gene_Mutation";

	// text field arrays holding the N probabilities
	double[] VD_N_Fields = null;

	double[] DJ_N_Fields = null;

	// text field arrays holding the Exonuclease probabilities
	// for Mean and Std Deviation for each exo region (start and end of genes)
	double[] V_end_exo_mean_Fields = null;

	double[] V_end_exo_stdDev_Fields = null;

	double[] D_start_exo_mean_Fields = null;

	double[] D_start_exo_stdDev_Fields = null;

	double[] D_end_exo_mean_Fields = null;

	double[] D_end_exo_stdDev_Fields = null;

	double[] J_start_exo_mean_Fields = null;

	double[] J_start_exo_stdDev_Fields = null;

	// text field arrays holding the P probabilities
	double[] V_end_P_Fields = null;

	double[] D_start_P_Fields = null;

	double[] D_end_P_Fields = null;

	double[] J_start_P_Fields = null;

	// text field holding the probability of a mutation in any nucleotide
	// irrelevant of Gene and Nucleotide position
	double Gene_Mutation_Field = 0.00;

	//************** Methods ******************//

	// create an array holding the exo nucl. probabilites based on a normal
	// distribution
	// the result is a normalised array ( the sum of the array equals ca. 1.0 )
	public double[] getExoProbArray(double mean, double stdDev,
			double minProbLimit, boolean reverseArray) {
		ArrayList temp_result = new ArrayList(); // temporary result
		double[] result_normalized;

		Double temp_D;
		double prob;

		for (int i = 0;; i++) {

			// calculate the probability of "i" exo nucl. removals
			prob = getNormalDistProb(i, mean, stdDev);
			temp_D = new Double(prob);

			if (prob < minProbLimit && i > mean) // test whether we are below
												 // min probability limit
				break; // we are finished
			else
				temp_result.add(temp_D);
		}

		// create an array holding the double values
		result_normalized = new double[temp_result.size()];

		//get the sum of the probabilities, so the results can be normalized
		double sum = 0.0;
		for (int i = 0; i < temp_result.size(); i++) {
			temp_D = (Double) temp_result.get(i);
			sum += temp_D.doubleValue();
		}

		// System.out.println("sum of probabilities = " + sum);

		// if the exo nucl. activity is at the end of a gene (3 prime), it is
		// easier to reverse the array
		if (reverseArray) {
			int array_index = 0;
			for (int i = (result_normalized.length - 1); i >= 0; i--) {
				temp_D = (Double) temp_result.get(i);
				result_normalized[array_index] = temp_D.doubleValue() / sum; // normalize
				++array_index;
			}
		} else // 5 prime
		{
			for (int i = 0; i < result_normalized.length; i++) {
				temp_D = (Double) temp_result.get(i);
				result_normalized[i] = temp_D.doubleValue() / sum; // normalize
			}
		}

		sum = 0.0;
		for (int i = 0; i < result_normalized.length; i++) {
			sum += result_normalized[i];
		}

		//	System.out.println("sum of normalized array = " + sum);

		return result_normalized;
	}//--getExoProbArray()

	private double getNormalDistProb(double x, double mean, double stdDev) {
		double result, temp1, temp2;

		temp1 = 1 / (Math.sqrt(2 * Math.PI * Math.pow(stdDev, 2)));
		temp2 = (Math
				.pow(Math.E, (-0.5 * (Math.pow(((x - mean) / stdDev), 2)))));
		result = temp1 * temp2;

		return result;
	}//--getNormalDistProb()

	public static ProbabilityHolder createFromString(String string) {
		String END = "end";

		ProbabilityHolder result = new ProbabilityHolder();

		StringTokenizer st = new StringTokenizer(string);

		try {
			String token, name_token;
			token = st.nextToken();

			if (!token.equalsIgnoreCase("Probabilities"))
				throw new Error();

			while (st.hasMoreTokens()) {
				token = st.nextToken();

				if (token.equals(VD_N)) {
					// System.out.println("parsing double array for VD_N_53");
					result.VD_N_Fields = parseDoubleArray(st);
				} else if (token.equals(DJ_N))
					result.DJ_N_Fields = parseDoubleArray(st);
				else if (token.equals(V_END_EXO_MEAN))
					result.V_end_exo_mean_Fields = parseDoubleArray(st);
				else if (token.equals(V_END_EXO_STDDEV))
					result.V_end_exo_stdDev_Fields = parseDoubleArray(st);
				else if (token.equals(D_START_EXO_MEAN))
					result.D_start_exo_mean_Fields = parseDoubleArray(st);
				else if (token.equals(D_START_EXO_STDDEV))
					result.D_start_exo_stdDev_Fields = parseDoubleArray(st);
				else if (token.equals(D_END_EXO_MEAN))
					result.D_end_exo_mean_Fields = parseDoubleArray(st);
				else if (token.equals(D_END_EXO_STDDEV))
					result.D_end_exo_stdDev_Fields = parseDoubleArray(st);
				else if (token.equals(J_START_EXO_MEAN))
					result.J_start_exo_mean_Fields = parseDoubleArray(st);
				else if (token.equals(J_START_EXO_STDDEV))
					result.J_start_exo_stdDev_Fields = parseDoubleArray(st);
				else if (token.equals(V_END_P))
					result.V_end_P_Fields = parseDoubleArray(st);
				else if (token.equals(D_START_P))
					result.D_start_P_Fields = parseDoubleArray(st);
				else if (token.equals(D_END_P))
					result.D_end_P_Fields = parseDoubleArray(st);
				else if (token.equals(J_START_P))
					result.J_start_P_Fields = parseDoubleArray(st);
				else if (token.equals(GENE_MUTATION))
					result.Gene_Mutation_Field = parseDouble(st);
			}//--while

			if (result.VD_N_Fields == null)
				System.out.println(" vd n 53 = null");
			else
				System.out.println(" vd n 53 not equal null");

		} catch (Exception e) {
			System.out.println("parsing exception");
		}

		return result;
	}//--createFromString

	/////////////////////////////////////////////////////////////

	private static double parseDouble(StringTokenizer st) {
		double d = Double.parseDouble(st.nextToken());
		System.out.println("parse double = " + d);
		return d;
	}//--parseDouble()

	/////////////////////////////////////////////////////////////

	private static double[] parseDoubleArray(StringTokenizer st) {
		ArrayList temp_result = new ArrayList();
		double[] result = null;
		Double temp_D;
		double temp_d;

		String token = st.nextToken();
		if (token.equalsIgnoreCase("Array")) {

			do {
				token = st.nextToken(); // test if this token is "end", else
										// skip it
				if (token.equals("end"))
					break;
				token = st.nextToken(); // this token should be a double
				temp_result.add(new Double(Double.parseDouble(token)));
			} while (st.hasMoreTokens());
		} else
			return null; // invalid parsing, throw an error in the future

		// create an array to return from the arraylist
		result = new double[temp_result.size()];

		for (int i = 0; i < result.length; i++) {
			temp_D = (Double) temp_result.get(i);
			temp_d = temp_D.doubleValue();
			result[i] = temp_d;
		}

		return result;
	}

	public String toString() {
		String N = "\n";
		String END = "end";
		String result = new String("Probabilities" + N + VD_N
				+ doubleArrayToString(VD_N_Fields) + END + N + DJ_N
				+ doubleArrayToString(DJ_N_Fields) + END + N +

				// text field arrays holding the Exonuclease probabilities
				// for Mean and Std Deviation for each exo region (start and end
				// of genes)
				V_END_EXO_MEAN + doubleArrayToString(V_end_exo_mean_Fields)
				+ END + N + V_END_EXO_STDDEV
				+ doubleArrayToString(V_end_exo_stdDev_Fields) + END + N
				+ D_START_EXO_MEAN
				+ doubleArrayToString(D_start_exo_mean_Fields) + END + N
				+ D_START_EXO_STDDEV
				+ doubleArrayToString(D_start_exo_stdDev_Fields) + END + N
				+ D_END_EXO_MEAN + doubleArrayToString(D_end_exo_mean_Fields)
				+ END + N + D_END_EXO_STDDEV
				+ doubleArrayToString(D_end_exo_stdDev_Fields) + END + N
				+ J_START_EXO_MEAN
				+ doubleArrayToString(J_start_exo_mean_Fields) + END + N
				+ J_START_EXO_STDDEV
				+ doubleArrayToString(J_start_exo_stdDev_Fields) + END + N +

				// text field arrays holding the P probabilities
				V_END_P + doubleArrayToString(V_end_P_Fields) + END + N
				+ D_START_P + doubleArrayToString(D_start_P_Fields) + END + N
				+ D_END_P + doubleArrayToString(D_end_P_Fields) + END + N
				+ J_START_P + doubleArrayToString(J_start_P_Fields) + END + N +

				// Gene Mutation Field
				GENE_MUTATION + "  " + Gene_Mutation_Field + "  " + END + N);

		return result;

	}//--toString()

	//////////////////////////////////////////////////////////////////////

	private String doubleArrayToString(double[] array) {
		String result = "  Array "; // initial spacing
		String temp;
		for (int i = 0; i < array.length; i++) {
			temp = "[" + i + "]  " + array[i] + "  ";
			result += temp;
		}

		return result;
	}

}//--ProbabilityHolder
